import {
  dsDNA,
  oligo,
  plasmid,
  polynucleotide,
  cleanup,
  resolveToSeq,
  revcomp,
  gccontent,
  basebalance,
  maxrepeat,
  isPalindromic
} from '../core/seq.js';

globalThis.dsDNA_JSON = (seq) => JSON.stringify(dsDNA(seq));
globalThis.oligo_JSON = (seq) => JSON.stringify(oligo(seq));
globalThis.plasmid_JSON = (seq) => JSON.stringify(plasmid(seq));
globalThis.polynucleotide_JSON = (...args) => JSON.stringify(polynucleotide(...args));

// Expose analysis tools directly
globalThis.cleanup = cleanup;
globalThis.resolveToSeq = resolveToSeq;
globalThis.revcomp = revcomp;
globalThis.gccontent = gccontent;
globalThis.basebalance = basebalance;
globalThis.maxrepeat = maxrepeat;
globalThis.isPalindromic = isPalindromic;
