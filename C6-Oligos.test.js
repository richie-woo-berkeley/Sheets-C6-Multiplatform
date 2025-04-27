import { describe, it, expect } from 'vitest';
import {
  isValidOligo,
  calculateTm,
  hasSelfComplementarity,
  designOligo,
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
  it('validates oligo sequences', () => {
    expect(isValidOligo('ATCG')).toBe(true);
    expect(isValidOligo('ATUG')).toBe(false);
  });

  it('calculates Tm correctly', () => {
    expect(calculateTm('ATGC')).toBe(2*2 + 4*2); // 2 A/T and 2 G/C
  });

  it('detects self-complementarity', () => {
    expect(hasSelfComplementarity('ATGCAA')).toBe(false);
  });

  it('designs oligos correctly', () => {
    expect(designOligo('ATGC', 'GAATTC')).toBe('GAATTCATGC');
  });

  it('scores annealing correctly', () => {
    const {score} = scoreanneal('ATGCATGTAAGTAATTTTAC', 'ATGTAAGTAATTTTAC');
    expect(score).toBeGreaterThan(10);
  });

  it('finds annealing positions', () => {
    const {position} = findanneal('ATGCATGTAAGTAATTTTACAGCTGGATTGCTATTACTTGTAATAGCATTTGGCGGAACATAA', 'ATGCATGTAAGTAATTTTAC');
    expect(position).toBeGreaterThanOrEqual(0);
  });

  it('finds partial complements with PCA', () => {
    expect(pca('CCATAGAATTCATGCCGTCTGGTCGGCCACTGGGGCACCACCCCGGGCCTGAACTTCCTTC', 'GATCAGGCGGTTGATGTGGGCGAGAAGGAAGTTCAGGCCCGGG')).toBeGreaterThan(0);
  });

  it('finds longest common annealing (LCA)', () => {
    expect(lca('CCATAGAATTCATGCCGTCTGGTCGGCCACTGGGGCACCACCCCGGGCCTGAACTTCCTTC', 'GATCAGGCGGTTGATGTGGGCGAGAAGGAAGTTCAGGCCCGGG')).toBeGreaterThan(0);
  });

  it('designs bglbrick oligos', () => {
    expect(bglbrick('ATGC')).toContain('AGATCT');
  });

  it('designs biobrick oligos', () => {
    expect(biobrick('ATGC')).toContain('GAATTCGCGGCCGCTTCTAGAG');
  });

  it('designs moclo oligos', () => {
    expect(moclo('ATGC')).toContain('GAAGAC');
  });

  it('designs genejoin oligos', () => {
    expect(genejoin('ATGC')).toContain('GGTCTC');
  });

  it('designs rbslib oligos', () => {
    expect(rbslib('ATGC')).toContain('AAGGAG');
  });
});
