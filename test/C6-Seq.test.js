import { describe, it, expect } from 'vitest';
import { cleanup, resolveToSeq, revcomp, isPalindromic, gccontent, basebalance, maxrepeat, translate, Polynucleotide, polynucleotide, comparePolynucleotides, plasmid, oligo, dsDNA, resolveToPoly } from 'src/C6-Seq.js';

describe('C6-Seq Utilities', () => {
  
  it('cleans up messy DNA input', () => {
    const messy = "4201 ACAACCCCAA GGACCGCGCC\nGAGAACCTCA TGATCGTGG";
    const cleaned = cleanup(messy);
    expect(cleaned).toBe('ACAACCCCAAGGACCGCGCCGAGAACCTCATGATCGTGG');
  });

  it('resolves simple valid sequences', () => {
    expect(resolveToSeq('atcg')).toBe('ATCG');
    expect(() => resolveToSeq('badseq!')).toThrow();
  });

  it('calculates correct reverse complement', () => {
    expect(revcomp('ATGC')).toBe('GCAT');
  });

  it('detects palindromic sequences', () => {
    expect(isPalindromic('AATT')).toBe(true);
    expect(isPalindromic('CGAT')).toBe(false);
    expect(() => isPalindromic('CGAN')).toThrow();
  });

  it('calculates GC content', () => {
    expect(gccontent('ACGT')).toBeCloseTo(0.5, 5);
    expect(gccontent('GGCC')).toBe(1);
  });

  it('calculates base balance', () => {
    expect(basebalance('ACGT')).toBeCloseTo(1.0, 5);
  });

  it('finds max repeat length', () => {
    expect(maxrepeat('AAATTTTGGG')).toBe(4); // 4 Ts
  });

  it('translates simple ORFs', () => {
    expect(translate('ATGTTT')).toBe('MF');
    expect(() => translate('ATGN')).toThrow();
  });

});

describe('Polynucleotide Objects', () => {

  it('creates a plasmid correctly', () => {
    const p = plasmid('ACAACCCCAAGGACCGGATCCGAGAACCTCATGATCGTGG');
    expect(p.isCircular).toBe(true);
    expect(p.isDoubleStranded).toBe(true);
  });

  it('creates an oligo correctly', () => {
    const o = oligo('ACAACCCCAAGGACCGGATCCGAGAACCTCATGATCGTGG');
    expect(o.isCircular).toBe(false);
    expect(o.isDoubleStranded).toBe(false);
  });

  it('creates a dsDNA fragment correctly', () => {
    const d = dsDNA('ACAACCCCAAGGACCGGATCCGAGAACCTCATGATCGTGG');
    expect(d.isCircular).toBe(false);
    expect(d.isDoubleStranded).toBe(true);
  });

  it('creates a sticky end fragment correctly', () => {
    const poly = polynucleotide('ACAACCCCAAGGACCGGATCCGAGAACCTCATGATCGTGG', 'AATT', 'GATC', true, false, false);
    expect(poly.ext5).toBe('AATT');
    expect(poly.ext3).toBe('GATC');
  });

  it('resolves from JSON or DNA string correctly', () => {
    const poly = resolveToPoly('ACAACCCCAAGGACCGGATCCGAGAACCTCATGATCGTGG');
    expect(poly.sequence).toBe('ACAACCCCAAGGACCGGATCCGAGAACCTCATGATCGTGG');
    expect(poly.isDoubleStranded).toBe(true);
  });

});
describe('comparePolynucleotides', () => {
  
  it('returns true for identical polynucleotides', () => {
    const poly1 = dsDNA('ATGCATGC');
    const poly2 = dsDNA('ATGCATGC');
    expect(comparePolynucleotides(poly1, poly2)).toBe(true);
  });

  it('returns false for different sequences', () => {
    const poly1 = dsDNA('ATGCATGC');
    const poly2 = dsDNA('ATGCATGA');
    expect(comparePolynucleotides(poly1, poly2)).toBe(false);
  });

  it('returns false for different 5\' overhangs', () => {
    const poly1 = polynucleotide('ATGCATGC', 'AATT', null, true, false, false);
    const poly2 = polynucleotide('ATGCATGC', 'TTAA', null, true, false, false);
    expect(comparePolynucleotides(poly1, poly2)).toBe(false);
  });

  it('returns false for different 3\' overhangs', () => {
    const poly1 = polynucleotide('ATGCATGC', null, 'GATC', true, false, false);
    const poly2 = polynucleotide('ATGCATGC', null, 'CTAG', true, false, false);
    expect(comparePolynucleotides(poly1, poly2)).toBe(false);
  });

  it('returns true when overhangs match but are empty', () => {
    const poly1 = polynucleotide('ATGCATGC', '', '', true, false, false);
    const poly2 = polynucleotide('ATGCATGC', '', '', true, false, false);
    expect(comparePolynucleotides(poly1, poly2)).toBe(true);
  });

});