# C6-Multiplatform

Cross-platform toolkit for synthetic biology and molecular biology operations.

[![npm version](https://img.shields.io/npm/v/c6-sim.svg?style=flat)](https://www.npmjs.com/package/c6-sim)
[![jsDelivr hits](https://data.jsdelivr.com/v1/package/npm/c6-sim/badge)](https://www.jsdelivr.com/package/npm/c6-sim)

---

## Overview

**C6** provides a suite of tools for parsing, simulating, and manipulating construction files (CFs) used in synthetic biology workflows.  
It includes functionality for:

- DNA sequence parsing and cleanup
- Molecular cloning simulations (PCR, Gibson assembly, Golden Gate)
- Restriction digest simulation
- Construction file parsing (text and spreadsheet-friendly formats)
- Sequence annotation and feature detection
- Oligonucleotide design and analysis
- Codon optimization and silent site removal

---

## Installation

Install via npm:

```bash
npm install c6-sim
```

Or load directly from CDN:

```html
<script src="https://cdn.jsdelivr.net/npm/c6-sim/dist/c6-sim.min.js"></script>
```

---

## Usage

After loading, access all functions via the global `C6` object (browser) or import (Node.js/ESM):

```javascript
// Example: Parse and simulate a construction file
const cf = C6.parseCF("PCR FwdPrimer RevPrimer on Template, Product");
cf.sequences = {
  "FwdPrimer": "ATGCGT",
  "RevPrimer": "TACGCA",
  "Template": "ATGCGTACGCATAGC"
};
const result = C6.simCF(cf);
console.log(result);
```

---

## Repository

[GitHub Repository](https://github.com/YOUR_USERNAME/C6-Multiplatform)

---

## License

MIT License © University of California, Berkeley
