"use strict";

/* ---------- FASTA parsing ---------- */
function parseFasta(text){
  const records = [];
  let header = null, seq = [];
  for (const raw of text.split(/\r?\n/)){
    const line = raw.trim();
    if (!line) continue;
    if (line[0] === ">"){
      if (header !== null) records.push({ header, seq: seq.join("") });
      header = line.slice(1).trim();
      seq = [];
    } else {
      seq.push(line.replace(/[^A-Za-z*]/g, "").toUpperCase());
    }
  }
  if (header !== null) records.push({ header, seq: seq.join("") });
  return records.filter(r => r.seq.length > 0);
}

/* match a FASTA header to one of the known proteins by alias */
function matchProtein(header){
  const h = header.toLowerCase();
  for (const [key, meta] of Object.entries(PROTEINS)){
    if (meta.aliases.some(a => h.includes(a))) return key;
  }
  return null;
}

/* ---------- sequence identity (banded Needleman–Wunsch, capped) ---------- */
/* Real global alignment identity, normalised over the aligned reference length.
   Sequences are capped so it stays fast in the browser. */
function identity(query, ref){
  const CAP = 500;
  const a = query.slice(0, CAP), b = ref.slice(0, CAP);
  const n = a.length, m = b.length;
  if (!n || !m) return 0;
  const GAP = -1, MATCH = 1, MIS = -1;
  let prev = new Int32Array(m + 1);
  let curr = new Int32Array(m + 1);
  // matches counter table via traceback would be heavy; approximate identity as
  // (score + m) / (2m) is unstable, so we do a lightweight matched-count DP:
  let prevMatch = new Int32Array(m + 1);
  let currMatch = new Int32Array(m + 1);
  for (let j = 0; j <= m; j++){ prev[j] = j * GAP; prevMatch[j] = 0; }
  for (let i = 1; i <= n; i++){
    curr[0] = i * GAP; currMatch[0] = 0;
    for (let j = 1; j <= m; j++){
      const eq = a[i-1] === b[j-1];
      const diag = prev[j-1] + (eq ? MATCH : MIS);
      const up   = prev[j] + GAP;
      const left = curr[j-1] + GAP;
      let best = diag, bm = prevMatch[j-1] + (eq ? 1 : 0);
      if (up > best){ best = up; bm = prevMatch[j]; }
      if (left > best){ best = left; bm = currMatch[j-1]; }
      curr[j] = best; currMatch[j] = bm;
    }
    [prev, curr] = [curr, prev];
    [prevMatch, currMatch] = [currMatch, prevMatch];
  }
  const matched = prevMatch[m];
  return Math.max(0, Math.min(1, matched / Math.min(n, m)));
}

/* ---------- scoring ---------- */
function scoreFromIdentities(idMap){
  const pairs = [];
  for (const p of PAIRS){
    const li = idMap[p.ligand], ri = idMap[p.receptor];
    if (li == null || ri == null) continue;              // pair not fully provided
    const limiting = Math.min(li, ri);                    // compatibility floor logic
    const cls = classifyIdentity(limiting);
    pairs.push({
      ligand: PROTEINS[p.ligand].label, receptor: PROTEINS[p.receptor].label,
      ligId: li, recId: ri, score: limiting, cls,
      limitedBy: li <= ri ? PROTEINS[p.ligand].label : PROTEINS[p.receptor].label,
    });
  }
  return pairs;
}

/* ---------- rendering ---------- */
const $ = id => document.getElementById(id);

function renderResults(species, pairs){
  $("emptyState").hidden = true;
  $("results").hidden = false;
  $("resultSpecies").textContent = species ? species : "";

  const order = { never:0, short:1, long:2 };
  const tracks = $("tracks");
  tracks.innerHTML = "";

  pairs.forEach((p, i) => {
    const row = document.createElement("div");
    row.className = "track";
    const flagged = p.cls === "never";
    row.innerHTML = `
      <div class="track-name">
        <span class="lig">${p.ligand}</span><span class="sep">·</span><span class="rec">${p.receptor}</span>
        ${flagged ? `<span class="track-bottleneck">bottleneck — limited by ${p.limitedBy}</span>` : ``}
      </div>
      <div class="lane">
        <span class="lane-tick t1"></span><span class="lane-tick t2"></span>
        <span class="lane-marker m-${p.cls}" style="left:0"></span>
      </div>
      <div class="track-id">${p.score.toFixed(2)}</div>`;
    tracks.appendChild(row);
    // animate marker to its ordinal position
    const marker = row.querySelector(".lane-marker");
    requestAnimationFrame(() => {
      setTimeout(() => { marker.style.left = (CLASS_POS[p.cls] * 100) + "%"; }, 40 + i * 70);
    });
  });

  // verdict: the limiting pair sets the maintenance ceiling
  const worst = pairs.reduce((w, p) => order[p.cls] < order[w.cls] ? p : w, pairs[0]);
  const flags = pairs.filter(p => p.cls === "never");
  const verdict = $("verdict");
  const ceiling = CLASS_LABEL[worst.cls];
  let line = `Predicted maintenance ceiling: <strong>${ceiling}</strong>. `;
  if (flags.length){
    line += `${flags.length} pair${flags.length>1?"s":""} below the compatibility floor — the rescue-cocktail target${flags.length>1?"s":""}.`;
  } else {
    line += `No pair falls below the floor.`;
  }
  verdict.innerHTML = `<p class="verdict-line">${line}</p>
    <div class="verdict-pairs">${
      pairs.map(p => `<span class="vchip ${p.cls==='never'?'flag':''}">${p.ligand}·${p.receptor} <b>${CLASS_LABEL[p.cls]}</b></span>`).join("")
    }</div>`;
}

/* ---------- controllers ---------- */
function runFromFasta(){
  const species = $("speciesName").value.trim();
  const text = $("fastaInput").value;
  const records = parseFasta(text);
  if (!records.length){
    setDetected("No FASTA records found. Paste sequences or upload a file.", true);
    return;
  }
  const idMap = {};
  const found = [];
  for (const r of records){
    const key = matchProtein(r.header);
    if (key && idMap[key] == null){
      idMap[key] = identity(r.seq, REFERENCE_SEQ[key]);
      found.push(PROTEINS[key].label);
    }
  }
  if (!found.length){
    setDetected("No records matched a known niche factor. Name headers e.g. >WNT3A, >FZD8.", true);
    return;
  }
  const pairs = scoreFromIdentities(idMap);
  if (!pairs.length){
    setDetected(`Matched ${found.join(", ")}, but no complete ligand·receptor pair. Provide both partners of a pair.`, true);
    return;
  }
  setDetected(`Matched ${found.length}: ${found.join(", ")}.`);
  renderResults(species, pairs);
}

function runZebrafishExample(){
  $("speciesName").value = "Danio rerio (zebrafish)";
  $("fastaInput").value = "";
  const idMap = {};
  for (const k of Object.keys(ZEBRAFISH_X1)) idMap[k] = ZEBRAFISH_X1[k];
  setDetected("Loaded real precomputed X1 identities for zebrafish (from the training set).");
  renderResults("Danio rerio", scoreFromIdentities(idMap));
}

function setDetected(msg, warn){
  const el = $("detected");
  el.hidden = false;
  el.textContent = msg;
  el.classList.toggle("warn", !!warn);
}

/* ---------- wire up ---------- */
$("predictBtn").addEventListener("click", runFromFasta);
$("exampleBtn").addEventListener("click", runZebrafishExample);
$("fastaFile").addEventListener("change", e => {
  const file = e.target.files[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = ev => {
    $("fastaInput").value = ev.target.result;
    const guess = file.name.replace(/\.(fasta|fa|faa|txt)$/i, "").replace(/[_-]+/g, " ");
    if (!$("speciesName").value) $("speciesName").value = guess;
    setDetected(`Loaded ${file.name}.`);
  };
  reader.readAsText(file);
});

// point the Source link at wherever this is hosted, if on GitHub
(function(){
  const link = document.getElementById("repoLink");
  if (location.hostname.endsWith("github.io")){
    const user = location.hostname.split(".")[0];
    const repo = location.pathname.split("/").filter(Boolean)[0] || "evocult-platform";
    link.href = `https://github.com/${user}/${repo}`;
  }
})();
