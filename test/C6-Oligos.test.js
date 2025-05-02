import { describe, it, expect } from 'vitest';
import {
  scoreanneal,
  findanneal,
  pca,
  lca,
  bglbrick,
  biobrick,
  moclo,
  genejoin,
  rbslib
} from 'src/C6-Oligos.js';


describe('C6-Oligos Utilities', () => {
  it('scores annealing correctly', () => {
    const score = scoreanneal('ATGCATGTAAGTAATTTTAC');
    expect(score).toBeGreaterThan(0);
    expect(score).toBeLessThanOrEqual(1);
  });

  it('finds annealing positions', () => {
    const best = findanneal('ATGCATGTAAGTAATTTTACAGCTGGATTGCTATTACTTGTAATAGCATTTGGCGGAACATAA', false, false);
    expect(best.length).toBeGreaterThanOrEqual(18);
    expect(best.length).toBeLessThanOrEqual(25);
  });

  it('finds partial complements with PCA', () => {
    const longSeq = 'ATGGCGTCTGGTCGACGTCGACGTCGACGTCGACGTCGACGTCGACGTCGAC'.repeat(2);
    const oligos = pca(longSeq);
    expect(Array.isArray(oligos)).toBe(true);
    expect(oligos.length).toBeGreaterThan(0);
  });

  it('finds longest common annealing (LCA)', () => {
    const oligos = lca('CCATAGAATTCATGCCGTCTGGTCGGCCACTGGGGCACCACCCCGGGCCTGAACTTCCTTC');
    expect(Array.isArray(oligos)).toBe(true);
    expect(oligos.length).toBeGreaterThan(0);
  });

  it('designs bglbrick oligos', () => {
    expect(bglbrick('ATGC', 'F')).toContain('AGATCT');
  });

  it('designs biobrick oligos', () => {
    expect(biobrick('ATGC', true, 'F')).toContain('GAATTCgcggccgctTCTAG');
  });

  it('designs moclo oligos', () => {
    expect(moclo('ATGC', 'C', 'F')).toContain('GGTCTC');
  });

  it('designs genejoin oligos', () => {
    expect(genejoin('ATGCATGC', 'TTAACCGG', 'F')).toMatch(/[ATGC]+/);
  });

  it('designs rbslib oligos', () => {
    expect(rbslib('ATGAAA', 'AGGAGGTTTA', 'F')).toMatch(/GGTCTCaTACT[a-z]+ATG/);
  });
});
