import { describe, it, expect } from 'vitest';
import { removeSites, oneAAoneCodon } from '../src/C6-Gene.js';
import { translate } from 'src/C6-Seq.js';

describe('C6-Gene Utilities', () => {
  it('translates a DNA sequence to protein correctly', () => {
    const dna = 'ATGTTCGGTCTCAACGGAGACCAGCAGGAATCTTAA';
    const protein = translate(dna);
    expect(protein).toBe('MFGLNGDQQES');
  });

  it('reverse-translates a protein to DNA correctly', () => {
    const protein = 'MFGLNGDQQES';
    const dna = oneAAoneCodon(protein);
    expect(dna).toBe('ATGTTTGGTTTAAATGGTGATCAACAAGAATCT');
    expect(translate(dna)).toBe(protein);
  });

  it('removes restriction sites from an ORF', () => {
    const orf = 'ATGTTCGGTCTCAACGGAGACCAGCAGGAATCTTAA';
    const cleaned = removeSites(orf);
    const cleanedProtein = translate(cleaned);
    expect(cleanedProtein).toBe('MFGLNGDQQES');
  });

  it('reverse-translates and removes sites for a protein sequence', () => {
    const protein = 'MDDASPRF*';
    const dna = oneAAoneCodon(protein);
    const cleaned = removeSites(dna);
    const translated = translate(cleaned);
    expect(translated).toBe('MDDASPRF');
  });
});
