import { resolveToSeq, gccontent, basebalance, maxrepeat, revcomp } from './seq.js';

export function scoreanneal(inseq) {
  const anneal = resolveToSeq(inseq);
  let score = 0;
  const maxPossibleScore = 5;

  if (['G', 'C'].includes(anneal[0])) score++;
  if (['G', 'C'].includes(anneal[anneal.length - 1])) score++;

  const gc = gccontent(anneal);
  if (gc >= 0.5 && gc <= 0.65) score++;

  const balance = basebalance(anneal);
  if (balance > 0.75) score++;

  const repeat = maxrepeat(anneal);
  if (repeat <= 3) score++;

  const lengthDiff = Math.abs(anneal.length - 20);
  score -= lengthDiff / 2;

  return Math.max(0, score / maxPossibleScore);
}

export function findanneal(inseq, lock5, lock3) {
  const seq = resolveToSeq(inseq);
  const minLength = 18;
  const maxLength = 25;
  let bestAnneal = "N/A";
  let bestScore = -1;

  if (lock5 && lock3) {
    throw new Error(`Cannot lock both ends of the template`);
  }

  if (lock5) {
    for (let end = minLength; end <= maxLength; end++) {
      const anneal = seq.substring(0, end);
      const score = scoreanneal(anneal);
      if (score > bestScore) {
        bestAnneal = anneal;
        bestScore = score;
      }
    }
    return bestAnneal;
  }

  if (lock3) {
    const end = seq.length;
    for (let start = end - maxLength; start < end - minLength; start++) {
      const anneal = seq.substring(start, end);
      const score = scoreanneal(anneal);
      if (score > bestScore) {
        bestAnneal = anneal;
        bestScore = score;
      }
    }
    return bestAnneal;
  }

  for (let start = 0; start < seq.length - minLength; start++) {
    for (let end = start + minLength; end <= Math.min(seq.length, start + maxLength); end++) {
      const anneal = seq.substring(start, end);
      const score = scoreanneal(anneal);
      if (score > bestScore) {
        bestAnneal = anneal;
        bestScore = score;
      }
    }
  }

  return bestAnneal;
}

// Return melting temperature
export function Tm(seq) {
  const s = resolveToSeq(seq);
  let tm = 0;
  const upper = s.toUpperCase();
  for (const base of upper) {
    if (base === 'G' || base === 'C') tm += 4;
    if (base === 'A' || base === 'T') tm += 2;
  }
  return tm;
}

// Return self-complementarity score
export function selfscore(seq) {
  const s = resolveToSeq(seq);
  const rev = revcomp(s);
  let score = 0;
  for (let i = 0; i < s.length; i++) {
    if (s[i] === rev[i]) score++;
  }
  return score;
}

// Design a forward primer
export function primerF(template, added5 = "") {
  const core = findanneal(template, true, false);
  return added5 + core;
}

// Design a reverse primer
export function primerR(template, added3 = "") {
  const core = findanneal(template, false, true);
  return revcomp(core + added3);
}

export function pca(synthon) {
  synthon = resolveToSeq(synthon);
  let chunks = Math.round(synthon.length / 25);
  if (chunks % 2 !== 0) chunks++;
  const chunksize = Math.round(synthon.length / chunks);

  const annealingIndices = [];
  for (let i = chunksize; i < synthon.length - chunksize; i += chunksize) {
    const seq = synthon.substr(i - 12, 24);
    const anneal = findanneal(seq, false, false);
    const startIndex = synthon.indexOf(anneal);
    const endIndex = startIndex + anneal.length;
    annealingIndices.push([startIndex, endIndex]);
  }

  const oligos = [];
  for (let i = 0; i < annealingIndices.length; i++) {
    if (i === 0) {
      oligos.push(synthon.substring(0, annealingIndices[1][1]));
    } else if (i === annealingIndices.length - 1) {
      oligos.push(revcomp(synthon.substring(annealingIndices[i][0])));
    } else if (i % 2 === 0) {
      oligos.push(synthon.substring(annealingIndices[i][0], annealingIndices[i + 1][1]));
    } else {
      oligos.push(revcomp(synthon.substring(annealingIndices[i][0], annealingIndices[i + 1][1])));
    }
  }

  return oligos;
}

export function lca(synthon) {
  synthon = resolveToSeq(synthon);
  const seqLen = synthon.length;
  const oligos = [];

  let mod = (seqLen - 25) % 50;
  let n = (seqLen - mod) / 50;
  if (mod > 25) n++;
  else if (mod < -25) n--;
  const chunkSize = Math.floor(seqLen / n);

  const generateOligos = (s) => {
    for (let i = 0; i < n; i++) {
      oligos.push(s.substring(i * chunkSize, i * chunkSize + chunkSize));
    }
  };

  generateOligos(synthon);
  generateOligos(revcomp(synthon));

  return oligos;
}

export function bglbrick(sequence, frgs) {
  const rORf = frgs[0].toUpperCase();
  sequence = resolveToSeq(sequence);

  if (rORf === 'F') {
    return "ccata" + "AGATCT" + findanneal(sequence, true, false);
  } else if (rORf === 'R') {
    return "catca" + "CTCGAGttaGGATCC" + revcomp(findanneal(sequence, false, true));
  } else if (rORf === 'S') {
    return "GAATTCatgAGATCT" + sequence + "GGATCCtaaCTCGAG";
  } else if (rORf === 'G') {
    return "ccataGAATTCatgAGATCT" + sequence + "GGATCCtaaCTCGAGtaacg";
  } else {
    throw new Error("Invalid value for 'frgs'. Must be 'F', 'R', 'S', or 'G'.");
  }
}

export function biobrick(sequence, isCDS, frgs) {
  const rORf = frgs[0].toUpperCase();
  sequence = resolveToSeq(sequence);

  if (rORf === 'F') {
    return isCDS
      ? "gacttGAATTCgcggccgctTCTAG" + findanneal(sequence, true, false)
      : "gacttGAATTCgcggccgctTCTAGAg" + findanneal(sequence, true, false);
  } else if (rORf === 'R') {
    return "catcaACTAGTa" + revcomp(findanneal(sequence, false, true));
  } else if (rORf === 'G') {
    return "ccataGAATTCgcggccgctTCTAG" + sequence + "tACTAGTagcggccgCTGCAGcatcg";
  } else if (rORf === 'S') {
    return isCDS
      ? "GAATTCgcggccgctTCTAG" + sequence + "tACTAGTagcggccgCTGCAG"
      : "GAATTCgcggccgctTCTAGAg" + sequence + "tACTAGTagcggccgCTGCAG";
  } else {
    throw new Error("Invalid value for 'frgs'. Must be 'F', 'R', 'S', or 'G'.");
  }
}

const stickyEnds = {
  UC: ['TACT', 'AAGC'],
  TP: ['GCTT', 'AGTA'],
  SP: ['AATG', 'ACCT'],
  U: ['TACT', 'CATT'],
  C: ['AGGT', 'AAGC'],
  T: ['GCTT', 'AGCG'],
  P: ['GGAG', 'AGTA']
};

export function moclo(sequence, partType, frgs) {
  const rORf = frgs[0].toUpperCase();
  sequence = resolveToSeq(sequence);
  const [sticky5, sticky3] = stickyEnds[partType] || [];

  if (!sticky5 || !sticky3) throw new Error("Invalid part type for MoClo.");

  if (rORf === 'F') {
    return "ccataGGTCTCa" + sticky5 + findanneal(sequence, true, false);
  } else if (rORf === 'R') {
    return "catcaGGTCTCt" + sticky3 + revcomp(findanneal(sequence, false, true));
  } else if (rORf === 'S') {
    return "GGTCTCt" + sticky5 + sequence + revcomp(sticky3) + "aGAGACC";
  } else if (rORf === 'G') {
    return "ccataGGTCTCt" + sticky5 + sequence + revcomp(sticky3) + "aGAGACCtaacg";
  } else {
    throw new Error("Invalid value for 'frgs'. Must be 'F', 'R', 'S', or 'G'.");
  }
}

export function genejoin(fivePrimeSeq, threePrimeSeq, ForR) {
  const seq5 = resolveToSeq(fivePrimeSeq);
  const seq3 = resolveToSeq(threePrimeSeq);
  const rORf = ForR[0].toUpperCase();

  const anneal5 = findanneal(seq5, false, true);
  const anneal3 = findanneal(seq3, true, false);
  const forOligo = anneal5 + anneal3;

  if (rORf === 'F') return forOligo;
  if (rORf === 'R') return revcomp(forOligo);

  throw new Error("Invalid value for 'ForR'. Must be 'F' or 'R'.");
}

export function rbslib(orf, utr, frg) {
  orf = resolveToSeq(orf);
  utr = resolveToSeq(utr);

  if (orf.length % 3 !== 0) throw new Error("Length of orf must be a multiple of 3");
  if (!["TAA", "TGA", "TAG"].includes(orf.slice(-3))) orf += "TAA";
  if (!["ATG", "GTG", "CTG", "TTG"].includes(orf.slice(0, 3))) {
    throw new Error("Start of CDS not a start codon");
  }
  if (!orf.startsWith("ATG")) orf = "A" + orf.substring(1);

  const rORf = frg[0].toUpperCase();
  const rbs = utr.slice(-17, -14) + "NVWGGRD" + utr.slice(-7);

  if (rORf === 'R') {
    return "catcaGGTCTCtAAGC" + revcomp(findanneal(orf, false, true));
  } else if (rORf === 'F') {
    return "ccataGGTCTCaTACT" + rbs.toLowerCase() + findanneal(orf, true, false);
  } else if (rORf === 'G') {
    return "ccataGGTCTCtTACT" + rbs.toLowerCase() + orf + "GCTTaGAGACCtgatg";
  } else if (rORf === 'S') {
    return "GGTCTCtTACT" + rbs.toLowerCase() + orf + "GCTTaGAGACC";
  }

  throw new Error("Invalid value for 'frg'. Must be 'F', 'R', 'G', or 'S'.");
}
