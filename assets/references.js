/* EvoCult reference panel ----------------------------------------------------
 * The eight core niche factors and their receptors, header-matching aliases,
 * the ligand-receptor pairs scored, and the compatibility floors.
 *
 * IMPORTANT: REFERENCE_SEQ entries below are ILLUSTRATIVE PLACEHOLDERS used only
 * so the browser prototype can compute a sequence-identity proxy out of the box.
 * In production, replace them with the curated human UniProt interface reference
 * panel (or move identity computation to the Python backend, which does this
 * against the real X1 pipeline). The zebrafish example does NOT use these — it
 * uses the real precomputed X1 identity values from the training set.
 * -------------------------------------------------------------------------- */

const PROTEINS = {
  // ligands
  WNT3A:  { role:"ligand",   label:"WNT3A",  aliases:["wnt3a","wnt-3a","wnt3"] },
  RSPO1:  { role:"ligand",   label:"RSPO1",  aliases:["rspo1","rspo-1","rspondin1","r-spondin1","r-spondin-1"] },
  EGF:    { role:"ligand",   label:"EGF",    aliases:["egf","epidermal growth factor"] },
  NOG:    { role:"ligand",   label:"Noggin", aliases:["nog","noggin"] },
  // receptors
  FZD8:   { role:"receptor", label:"FZD8",   aliases:["fzd8","frizzled8","frizzled-8","fzd-8"] },
  LGR5:   { role:"receptor", label:"LGR5",   aliases:["lgr5","lgr-5"] },
  EGFR:   { role:"receptor", label:"EGFR",   aliases:["egfr","egf receptor","erbb1","her1"] },
  BMPR1A: { role:"receptor", label:"BMPR1A", aliases:["bmpr1a","bmpr-1a","alk3"] },
  BMPR1B: { role:"receptor", label:"BMPR1B", aliases:["bmpr1b","bmpr-1b","alk6"] },
  LRP6:   { role:"receptor", label:"LRP6",   aliases:["lrp6","lrp-6"] },
};

// ligand · receptor pairs evaluated (from the zebrafish rescue design)
const PAIRS = [
  { ligand:"WNT3A", receptor:"FZD8"   },
  { ligand:"EGF",   receptor:"EGFR"   },
  { ligand:"RSPO1", receptor:"LGR5"   },
  { ligand:"NOG",   receptor:"BMPR1A" },
  { ligand:"NOG",   receptor:"BMPR1B" },
];

// compatibility floors on the limiting (min) pair identity — illustrative,
// calibrated from the training-set X1 distribution. Production uses the empirical
// floor of the validated winning feature instead of raw sequence identity.
const FLOORS = { never_below: 0.55, long_from: 0.75 };

function classifyIdentity(v){
  if (v < FLOORS.never_below) return "never";
  if (v >= FLOORS.long_from)  return "long";
  return "short";
}
const CLASS_LABEL = { never:"Never", short:"Short-term", long:"Long-term" };
// ordinal position 0..1 for the marker, placed mid-band
const CLASS_POS = { never:0.16, short:0.5, long:0.84 };

// real precomputed X1 interface sequence-identity values for zebrafish
// (Danio Rerio), taken directly from the EvoCult training table.
const ZEBRAFISH_X1 = {
  WNT3A:0.8958, EGF:0.4198, RSPO1:0.5364, NOG:0.6099,
  FZD8:0.6729,  EGFR:0.646, LGR5:0.5681,  BMPR1A:0.7752, BMPR1B:0.8313, LRP6:0.8158,
};

/* Illustrative placeholder references (see header note). Short canonical-region
 * fragments — enough for the proxy to run; NOT authoritative sequences. */
const REFERENCE_SEQ = {
  WNT3A:"MAPLGYFLLLCSLKQALGSYPIWWSLAVGPQYSSLGSQPILCASIPGLVPKQLRFCRNYVEIMPSVAEGIKIGIQECQHQFRGRRWNCTTVHDSLAIFGPVLDKATRESAFVHAIASAGVAFAVTRSCAEGTAAICGCSSRHQGSPGKGWKWGGCSEDIEFGGMVSREFADARENRPDARSAMNRHNNEAGRQAIASHMHLKCKCHGLSGSCEVKTCWWSQPDFRAIGDFLKDKYDSASEMVVEKHRESRGWVETLRPRYTYFKVPTERDLVYYEASPNFCEPNPETGSFGTRDRTCNVSSHGIDGCDLLCCGRGHNARAERRREKCRCVFHWCCYVSCQECTRVYDVHTCK",
  RSPO1:"MRLGLCVVALVLSWTHLTISSRGIKGKRQRRISAEGSQACAKGCELCSEVNGCLKCSPKLFILLERNDIRQVGVCLPSCPPGYFDARNPDMNKCIKCKIEHCEACFSHNFCTKCKEGLYLHKGRCYPACPEGSSAANGTMECSSPAQCEMSEWSPWGPCSKKQQLCGFRRGSEERTRRVLHAPVGDHAACSDTKETRRCTVRRVPCPEGQKRRKGGQGRRENANRNLARKESKEAGAGSRRRKGQQQQQQQGTVGPLTSAGPA",
  EGF:"NSDSECPLSHDGYCLHDGVCMYIEALDKYACNCVVGYIGERCQYRDLKWWELR",
  NOG:"MERCPSLGVTLYALVVVLGLRATPAGGQHYLHIRPAPSDNLPLVDLIEHPDPIFDPKEKDLNETLLRSLLGGHYDPGFMATSPPEDRPGGGGGAAGGAEDLAELDQLLRQRPSGAMPSEIKGLEFSEGLAQGKKQRLSKKLRRKLQMWLWSQTFCPVLYAWNDLGSRFWPRYVKVGSCFSKRSCSVPEGMVCKPSKSVHLTVLRWRCQRRGGQRCGWIPIQYPIISECKCSC",
  FZD8:"MEWGYLLEVTSLLAALALLQRSSGAAAASAKELACQEITVPLCKGIGYNYTYMPNQFNHDTQDEAGLEVHQFWPLVEIQCSPDLKFFLCSMYTPICLEDYKKPLPPCRSVCERAKAGCAPLMRQYGFAWPDRMRCDRLPEQGNPDTLCMDYNRTDLTTAAPSPPRRLPPPPPGEQPPSGSGHGRPPGARPPHRGGGRGGGGDAAAPPARGGGGGGKARPPGGGAAPCEPGCQCRAPMVSVSSERHPLYNRVKTGQIANCALPCHNPFFSQDERAFTVFWIGLWSVLCFVSTFATVSTFLIDMERFKYPERPIIFLSGCYFMVAVAHVAGFLLEDRAVCVERFSDDGYRTVAQGTKKEGCTILFMMLYFFGMASSIWWVILss",
  LGR5:"MDTSRLGVLLSLPVLLQLATGGSSPRSGVSCPQECTCQPTPMTVDCHELGLTAVPTNIPADTRHLNLQ",
  EGFR:"LEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYVVTDHGSCVRACGADSYEMEEDGVRKCKKCEGPCRKV",
  BMPR1A:"MPQLYIYIRLLGAYLFIISRVQGQNLDSMLHGTGMKSDSDQKKSENGVTLAPEDTLPFLKCYCSGHCPDDAINNTCITNGHCFAIIEEDDQGETTLASGCMKYEGSDFQCKDSPKAQLRRTIECCRTNLCNQYLQPTLPPVVIGPFFDGSIR",
  BMPR1B:"LLRSSPGKQTFRRTQERPVTPFPCLLLLLLPPTPGCGARIPTLPPPAPLLPHLLPTLLLLLPPPPLPVPPPPPADGADQVQKKNAIHTTFTCPTQSPQ",
  LRP6:"MGAVLRSLLACSFCVLLRAAPLLLYANRRDLRLVDATNGKENATIVVGGLEDAAAVDFVFSHGLIYWSDVSEEAIKRTEFNKTESVQNVVVSGLLSPDGLACDWLGEKLYWTDSETNRIEVANLNGTSRKVLFWQDLDQPRAIALDPAHGYMYWTDWGETPRIERAGMDGSTRKIIVDSDIYWPNGLTIDLEEQKLYWADAKLSFIHRANLDGSFRQKVVEGSLTHPFALTLSGDTLYWTDWQTRSIHACNKRTGGKRKEILSALYSPMDIQVLSQERQP",
};
