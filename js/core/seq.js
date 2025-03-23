export function cleanup(sequence) {
    if (typeof sequence !== 'string') sequence = sequence.toString();
    sequence = sequence.replace(/[\s\d\r]/g, "");
    const validBiopolymer = /^[ACGTRYSWKMBDHVNUacgtryswkmbdhvnu*]+$/;
    if (!validBiopolymer.test(sequence)) {
        throw new Error("Invalid characters in sequence.");
    }
    return sequence.toUpperCase();
}

const _regexDNA = /^[ACTGactgMRWSYKVHDBNXmrwsykvhdbnx-]+$/;
export function resolveToSeq(seq) {
    seq = seq.toString();
    if (_regexDNA.test(seq)) return seq.toUpperCase();
    throw new Error("Unrecognizable as sequence: " + seq);
}

export class Polynucleotide {
    constructor(sequence, ext5 = null, ext3 = null, isDoubleStranded = true, isRNA = false, isCircular = false, mod_ext5 = "", mod_ext3 = "") {
        this.sequence = sequence.toUpperCase();
        this.ext5 = ext5 ? ext5.toUpperCase() : ext5;
        this.ext3 = ext3 ? ext3.toUpperCase() : ext3;
        this.isDoubleStranded = isDoubleStranded;
        this.isRNA = isRNA;
        this.isCircular = isCircular;
        this.mod_ext5 = mod_ext5;
        this.mod_ext3 = mod_ext3;
    }

    toJSON() {
        return {
            sequence: this.sequence,
            ext5: this.ext5,
            ext3: this.ext3,
            isDoubleStranded: this.isDoubleStranded,
            isRNA: this.isRNA,
            isCircular: this.isCircular,
            mod_ext5: this.mod_ext5,
            mod_ext3: this.mod_ext3
        };
    }
}

export function _resolveToPoly(seqOrJSON) {
    try {
        return JSON.parse(seqOrJSON);
    } catch (err) {
        if (!_regexDNA.test(seqOrJSON)) {
            throw new Error("Cannot resolve " + seqOrJSON);
        }
        return new Polynucleotide(seqOrJSON, "", "", true, false, false, "hydroxyl", "hydroxyl");
    }
}

export function dsDNA(sequence) {
    return new Polynucleotide(sequence, "", "", true, false, false, "hydroxyl", "hydroxyl");
}
export function oligo(sequence) {
    return new Polynucleotide(sequence, null, null, false, false, false, "hydroxyl", "");
}
export function plasmid(sequence) {
    return new Polynucleotide(sequence, "", "", true, false, true, "", "");
}

export function revcomp(seq) {
    const map = {
        A: 'T', T: 'A', C: 'G', G: 'C',
        a: 't', t: 'a', c: 'g', g: 'c',
        B: 'V', D: 'H', H: 'D', K: 'M', N: 'N',
        R: 'Y', S: 'S', M: 'K', V: 'B', W: 'W', Y: 'R',
        b: 'v', d: 'h', h: 'd', k: 'm', n: 'n',
        r: 'y', s: 's', m: 'k', v: 'b', w: 'w', y: 'r'
    };
    return [...seq].reverse().map(c => map[c] || (() => { throw new Error(`Invalid base '${c}'`) })()).join('');
}

export function gccontent(seq) {
    const upper = seq.toUpperCase();
    const gc = (upper.match(/[GC]/g) || []).length;
    return gc / upper.length;
}

export function isPalindromic(seq) {
    return seq === revcomp(seq);
}

export function basebalance(seq) {
    const counts = { A: 0, C: 0, G: 0, T: 0 };
    const upper = seq.toUpperCase();
    for (const base of upper) if (counts[base] !== undefined) counts[base]++;
    let score = 1;
    for (const base in counts) {
        if (counts[base] === 0) return 0;
        score *= counts[base] / seq.length;
    }
    return 4 * Math.pow(score, 1/4);
}

export function maxrepeat(seq) {
    const upper = seq.toUpperCase();
    let max = 0, count = 0, last = '';
    for (const base of upper) {
        if (base === last) {
            count++;
            max = Math.max(max, count);
        } else {
            count = 1;
            last = base;
        }
    }
    return max;
}

export function polynucleotide(sequence, ext5, ext3, isDoubleStranded, isRNA, isCircular, mod_ext5, mod_ext3) {
    return new Polynucleotide(sequence, ext5, ext3, isDoubleStranded, isRNA, isCircular, mod_ext5, mod_ext3);
}
