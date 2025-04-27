import { annotateSequenceSmart, inferTranscriptionalUnits, inferExpressedProteins, findNonExpressedCDS } from 'src/C6-Annotator.js';

describe('C6-Annotator basic tests', () => {
  test('annotateSequenceSmart finds T7 promoter', async () => {
    const fakeFeatureDb = [
      { Name: 'T7 Promoter', Sequence: 'TAATACGACTCACTATAGGG', Type: 'promoter', Color: '#ff0000' }
    ];
    const sequence = 'GAAATTAATACGACTCACTATAGGGGAAT';
    const annotations = annotateSequenceSmart(sequence, fakeFeatureDb);

    expect(annotations.length).toBeGreaterThan(0);
    expect(annotations[0].label).toBe('T7 Promoter');
    expect(annotations[0].strand).toBe(1);
  });

  test('inferTranscriptionalUnits returns empty list for empty features', () => {
    const tus = inferTranscriptionalUnits([]);
    expect(tus).toEqual([]);
  });

  test('inferExpressedProteins returns empty array if no CDS', () => {
    const tus = [{ features: [{ type: 'promoter', label: 'p1' }] }];
    const proteins = inferExpressedProteins(tus);
    expect(proteins.length).toBe(0);
  });
  
  test('findNonExpressedCDS correctly identifies no non-expressed CDS in a multi-cds TU', () => {
    const allFeatures = [
      { type: 'cds', label: 'Protein1' },
      { type: 'cds', label: 'Protein2' }
    ];
    const expressedProteins = [
      { label: 'Protein1' },
      { label: 'Protein2' }
    ];
    const nonExpressed = findNonExpressedCDS(allFeatures, expressedProteins);
    expect(nonExpressed.length).toBe(0);  // No non-expressed CDS
  });

});

describe('C6-Annotator synthetic TU and expression tests', () => {
  test('TU inference detects simple promoter-cds-terminator', () => {
    const features = [
      { type: 'promoter', label: 'P1', start: 0, end: 10 },
      { type: 'cds', label: 'Gene1', start: 11, end: 50 },
      { type: 'terminator', label: 'T1', start: 51, end: 60 }
    ];
    const tus = inferTranscriptionalUnits(features);
    expect(tus.length).toBe(1);
    expect(tus[0].promoter.label).toBe('P1');
    expect(tus[0].terminator.label).toBe('T1');
    expect(tus[0].features.length).toBe(1);
    expect(tus[0].features[0].label).toBe('Gene1');
  });

  test('expressed proteins detected from TU', () => {
    const sequence = 'ATGCGT';
    const tus = [{ features: [{ type: 'cds', label: 'Gene1' }] }];
    const proteins = inferExpressedProteins(tus);
    expect(proteins.length).toBe(1);
    expect(proteins[0].label).toBe('Gene1');
  });

  test('finds multiple TUs and handles non-expressed CDS', () => {
    const features = [
      { type: 'promoter', label: 'P1', start: 0, end: 10 },
      { type: 'cds', label: 'Gene1', start: 11, end: 50 },
      { type: 'terminator', label: 'T1', start: 51, end: 60 },
      { type: 'cds', label: 'Gene2', start: 70, end: 100 }
    ];
    const tus = inferTranscriptionalUnits(features);
    const proteins = inferExpressedProteins(tus);
    const nonExpressed = findNonExpressedCDS(features, proteins);

    expect(tus.length).toBe(1);
    expect(proteins.length).toBe(1);
    expect(nonExpressed.length).toBe(1);
    expect(nonExpressed[0].label).toBe('Gene2');
  });
});