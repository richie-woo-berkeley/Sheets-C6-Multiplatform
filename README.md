# C6-Multiplatform

**Cross-platform toolkit for synthetic biology and molecular biology operations**

[![npm version](https://img.shields.io/npm/v/c6-sim.svg?style=flat)](https://www.npmjs.com/package/c6-sim)  
[![jsDelivr hits](https://data.jsdelivr.com/v1/package/npm/c6-sim/badge)](https://www.jsdelivr.com/package/npm/c6-sim)

---

## Overview

C6 is a JavaScript toolkit for simulating recombinant DNA workflows. It supports both browser-based and Node.js-based execution environments, enabling use in web applications, teaching tools, or backend design pipelines.

The central abstraction in C6 is the **Construction File (CF)** — a stepwise text specification describing a DNA assembly plan using PCR, restriction digest, Gibson assembly, Golden Gate cloning, transformation, and other procedures. C6 parses CFs into structured instructions and simulates the resulting chemical outcomes using modeled reactions.
The Construction File format is introduced in [Anderson et al., 2023](https://pubmed.ncbi.nlm.nih.gov/37956196/), which provides the conceptual foundation for this abstraction.

C6 models DNA molecules using a `Polynucleotide` abstraction that reflects the physical and chemical properties of DNA fragments used in recombinant DNA procedures. This includes parameters such as strandedness (single or double), topology (linear or circular), RNA vs DNA, and 5′/3′ overhangs or modifications. Helper constructors are provided for common cases (e.g., `dsDNA`, `oligo`, `plasmid`), but the `Polynucleotide` class itself generalizes any linear or circular molecule encountered in assembly workflows.

Beyond simulation, C6 includes design and verification tools:

- Oligo design utilities for standard assembly schemes (e.g., BioBrick, MoClo, BglBrick, RBSLib)
- Annotation and inference tools that map known features to sequences, group them into transcriptional units (TUs), and infer which CDS elements are likely expressed
- Sequence manipulation utilities for reverse complementation, melting temperature, self-complementarity checks, and more


C6 focuses on modeling the chemical logic of molecular biology protocols rather than abstract logic of synthetic circuits. It is useful for applications involving validation, verification, or generation of DNA designs in formats that correspond closely to wet-lab workflows and in vitro molecular biology techniques.

C6 is not a compositional design tool. It does not attempt to optimize or select parts for a genetic circuit, metabolic pathway, or regulatory network. Instead, it assumes a design has already been made and focuses on the molecular construction strategy that will realize it. It verifies that a Construction File is syntactically well-formed, chemically coherent, and capable of producing the intended DNA product in silico. It also provides limited analysis of potential gene expression from the resulting constructs. While C6 does not yet generate detailed wet-lab protocols or robotic instructions, it offers a foundation for that by providing precise modeling of DNA species and cloning operations.

---

## Installation

```bash
npm install c6-sim
```

Or load via CDN:

```html
<script src="https://cdn.jsdelivr.net/npm/c6-sim/dist/c6-sim.min.js"></script>
```

---
## Quickstart

You can use C6 to simulate recombinant DNA steps in Node.js or in the browser.

```javascript
import C6 from 'c6-sim';

const cfText = `
PCR Fwd Rev Template Product
oligo Fwd CAAGTGGGAACGCGTAATG
oligo Rev CGGTCACGGCACCACCATC
plasmid Template TTCAAGTGGGAACGCGTAATGAATTTTGAAGATGGTGGTGCCGTGACCG
`;

const cf = C6.parseCF(cfText);
const result = C6.simCF(cf);
console.log(result);
```

---

## Interactive Demo

You can try the toolkit in your browser using the hosted demo:

👉 [**C6 Interactive Demo**](https://ucb-bioe-anderson-lab.github.io/C6-Multiplatform)

This includes CF simulation, annotation, and expression inference in a no-install playground.

---

## Running Tests

This project uses [Vitest](https://vitest.dev):

```bash
npm test
```

Test coverage includes core simulation, sequence utilities, oligo design, and gene annotation.

---

## Repository

[GitHub – UCB BioE Anderson Lab](https://github.com/UCB-BioE-Anderson-Lab/C6-Multiplatform)

---

## License

MIT License © University of California, Berkeley
