// C6-Annotator.js - Improved DNA Autoannotation and Transcription Inference Library for Web

import { cleanup, revcomp, translate } from './C6-Seq.js';


// Pre-index features by k-mers for fast lookup
function buildFeatureIndex(features, k = 10) {
  console.log("🔧 Building feature index...");
  const index = new Map();

  features.forEach(feature => {
    const pattern = cleanup(feature.Sequence || '');
    if (pattern.length >= k) {
      for (let i = 0; i <= pattern.length - k; i++) {
        const kmer = pattern.slice(i, i + k);
        if (!index.has(kmer)) {
          index.set(kmer, []);
        }
        index.get(kmer).push(feature);
      }
    }
  });

  console.log(`✅ Feature index built with ${index.size} unique ${k}-mers.`);
  return index;
}

// Smart matching using full exact matching (no k-mer seeding)
function annotateSequenceSmart(sequence, featureDb = null) {
  //   console.log("🔍 Starting annotation...");
  sequence = cleanup(sequence);
  const detectedFeatures = [];

  const db = featureDb || featureDbGlobal;
  if (!db) throw new Error("No feature database loaded yet.");

  const seqVariants = [sequence, revcomp(sequence)];
  //   console.log("🧬 Scanning sequence and reverse complement with full exact matching...");

  seqVariants.forEach((seq, strandIndex) => {
    db.forEach(feature => {
      const pattern = cleanup(feature.Sequence || '');
      if (pattern.length < 10) return; // Ignore very short patterns

      let pos = seq.indexOf(pattern);
      while (pos !== -1) {
        detectedFeatures.push({
          start: pos,
          end: pos + pattern.length,
          strand: strandIndex === 0 ? +1 : -1,
          label: feature.Name,
          type: feature.Type,
          color: feature.Color
        });
        pos = seq.indexOf(pattern, pos + 1);
      }
    });
  });

  //   console.log(`🔎 Found ${detectedFeatures.length} matching features.`);
  return detectedFeatures.sort((a, b) => a.start - b.start);
}

function inferTranscriptionalUnits(features) {
  //   console.log("🧬 Starting new-style TU inference...");

  const tus = [];
  const featureList = features.slice().sort((a, b) => a.start - b.start);
  const openTUs = [];

  const allowedTypes = new Set([
    'promoter', 'cds', 'terminator', 'rbs', 'kozak', 'polyA_signal',
    'intron', 'exon', 'utr', 'recombination_site', 'operator', 'enhancer',
    'silencer', 'riboswitch', 'insulator'
  ]);

  for (const feature of featureList) {
    const type = feature.type.toLowerCase();
    if (!allowedTypes.has(type)) continue;

    if (type === 'promoter') {
      //   console.log(`🔵 Found promoter: ${feature.label} at ${feature.start}`);
      openTUs.push({
        promoter: feature,
        start: feature.end, // Start at end of promoter
        features: [],
        terminator: null,
        end: null
      });
    }

    // Add feature to all open TUs
    openTUs.forEach(tu => {
      if (feature.start >= tu.start) {
        tu.features.push(feature);
        // console.log(`➕ Assigned feature ${feature.label} (${feature.type}) to TU started by ${tu.promoter.label}`);
      }
    });

    if (type === 'terminator') {
      //   console.log(`🔴 Found terminator: ${feature.label} at ${feature.start}`);
      // Close all open TUs
      openTUs.forEach(tu => {
        tu.terminator = feature;
        tu.end = feature.start; // End before terminator starts
        tus.push(tu);
        // console.log(`✅ Closed TU from ${tu.start} to ${feature.start} (promoter: ${tu.promoter.label}, terminator: ${feature.label})`);
      });
      openTUs.length = 0; // Clear open TUs
    }
  }

  //   console.log(`✅ Finished TU inference: ${tus.length} transcriptional units.`);
  return tus;
}


// Infer expressed proteins from TUs
function inferExpressedProteins(sequence, tus) {
  const proteins = [];

  tus.forEach((tu, index) => {
    tu.features.forEach(feature => {
      if (feature.type.toLowerCase() === 'cds') {
        proteins.push({
          tuIndex: index + 1,
          label: feature.label
        });
      }
    });
  });

  return proteins;
}

// Find non-expressed CDS
function findNonExpressedCDS(allFeatures, expressedProteins) {
  const expressedLabels = new Set(expressedProteins.map(p => p.label));
  const nonExpressed = [];

  allFeatures.forEach(feature => {
    if (feature.type.toLowerCase() === 'cds' && !expressedLabels.has(feature.label)) {
      nonExpressed.push({ label: feature.label });
    }
  });

  return nonExpressed;
}

// Internal feature database
let featureDbGlobal = [];

// Load feature database automatically
(function initializeFeatureDatabase() {
  const defaultFeatureUrl = "https://raw.githubusercontent.com/UCB-BioE-Anderson-Lab/cloning-tutorials/main/sequences/Default_Features.txt";

  //   console.log("🌐 Fetching default features...");
  fetch(defaultFeatureUrl)
    .then(response => {
      //   console.log("📥 Feature file fetched, parsing...");
      return response.text();
    })
    .then(text => {
      const lines = text.split("\n").filter(line => line.trim().length > 0);
      featureDbGlobal = lines.map(line => {
        const [Name, Sequence, Type, Color, LabelColor, Forward, Reverse] = line.split(/\s+/);
        return { Name, Sequence, Type, Color };
      });
      //   console.log(`✅ C6-Annotator: Loaded ${featureDbGlobal.length} features.`);
    })
    .catch(err => {
      console.error("❌ Failed to load default features:", err);
    });
})();

export {
  annotateSequenceSmart,
  inferTranscriptionalUnits,
  inferExpressedProteins,
  findNonExpressedCDS,
  featureDbGlobal
};