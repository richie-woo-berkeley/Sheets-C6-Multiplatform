import { describe, it, expect } from 'vitest';
import { PCR, assemble, gibson, cutOnce, digest, parseCF, simCF } from 'src/C6-Sim.js';

describe('C6-Sim Tests', () => {
  
  it('PCR correctly predicts PCR product', () => {
    const forward = 'gacttGAATTCgcggccgctTCTAGAgTCCCTATCAGTGATAGAG';
    const reverse = 'catcaACTAGTaGTGCTCAGTATCTCTATCAC';
    const template = 'tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac';
    const result = PCR(forward, reverse, template);
    expect(result).toBe('GACTTGAATTCGCGGCCGCTTCTAGAGTCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGATACTGAGCACTACTAGTTGATG');
  });

  it('assemble correctly assembles fragments by Golden Gate', () => {
    const frag1 = 'ccaaaGGTCTCAGCTTTGATCGATTCAACCTACTTCCCCTTCATAATCGGTACTAGAGACCacgac';
    const frag2 = 'GGTCTCATACTCAAAATTTACTGACTGGACATGGTCACCACTTAAGTAAGCTTTGAGACC';
    const result = assemble(frag1, frag2, 'BsaI');
    expect(result).toBe('GCTTTGATCGATTCAACCTACTTCCCCTTCATAATCGGTACTCAAAATTTACTGACTGGACATGGTCACCACTTAAGTAA');
  });

  it('gibson correctly assembles fragments by homologous recombination', () => {
    const frag1 = 'TTGGGTGCACGAGTGGGTTACatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtATTGACGCCGGGCAAGAGCAACT';
    const frag2 = 'ATTGACGCCGGGCAAGAGCAACTcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacATGGGGGATCATGTAACTCg';
    const frag3 = 'ATGGGGGATCATGTAACTCgccttgatcgaaagtaaaagatgctgaagatcagTTGGGTGCACGAGTGGGTTAC';
    const result = gibson([frag1, frag2, frag3]);
    expect(result).toMatch(/ATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGT/); // Start matches expected sequence
    expect(result.length).toBe(398);
  });

  it('cutOnce correctly cuts linear DNA with BamHI', () => {
    const poly = {
      sequence: "ACAACCCCAAGGACCGGATCCGAGACCCTGCAGTGATCGTGG",
      ext5: "",
      ext3: "",
      isDoubleStranded: true,
      isRNA: false,
      isCircular: false,
      mod_ext5: "hydroxyl",
      mod_ext3: "hydroxyl"
    };
    const result = cutOnce(poly, "BamHI");
    expect(result.length).toBe(2);
    expect(result[0].ext3).toBe('GATC');
    expect(result[1].ext5).toBe('GATC');
  });

  it('digest correctly digests a plasmid to completion', () => {
    const plasmid = {
      sequence: "AAAAAGAATTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGATCCGGGGG",
      ext5: "",
      ext3: "",
      isDoubleStranded: true,
      isRNA: false,
      isCircular: true,
      mod_ext5: "",
      mod_ext3: ""
    };
    const result = digest(plasmid, 'EcoRI,BamHI', 1);
    expect(result.sequence).toBe('CTTTTTTTTTTTTTTTTTTTTTTTTTTTTG');
    expect(result.ext5).toBe('AATT');
    expect(result.ext3).toBe('GATC');
  });

  it('parseCF parses simple construction file into correct steps and sequences', () => {
    const input = [
      ["PCR", "P6libF", "P6libR", "on", "pTP1", "P6"],
      ["Assemble", "P6", "BsaI", "pP6"],
    ];
    const parsed = parseCF(input);
    expect(parsed.steps.length).toBe(2);
    expect(parsed.steps[0].operation).toBe('PCR');
    expect(parsed.steps[1].operation).toBe('Assemble');
  });

  it('simCF runs construction steps and outputs correct product names', () => {
    const input = {
      steps: [
        {
          operation: "PCR",
          forward_oligo: "fwd",
          reverse_oligo: "rev",
          template: "tmpl",
          output: "pcrProduct"
        }
      ],
      sequences: {
        fwd: "gacttGAATTCgcggccgctTCTAGAgTCCCTATCAGTGATAGAG",
        rev: "catcaACTAGTaGTGCTCAGTATCTCTATCAC",
        tmpl: "tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac"
      }
    };
    const output = simCF(input);
    expect(output[0][0]).toBe('pcrProduct');
  });

});
