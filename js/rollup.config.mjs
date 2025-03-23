import { nodeResolve } from '@rollup/plugin-node-resolve';
import terser from '@rollup/plugin-terser';

export default [
  {
    input: 'adapters/web.js',
    output: {
      file: 'dist/c6-web.js',
      format: 'iife',
      name: 'C6Tools'
    },
    plugins: [nodeResolve(), terser()]
  },
  {
    input: 'adapters/sheets.js',
    output: {
      file: 'dist/c6-sheets.gs',
      format: 'iife'
    },
    plugins: [nodeResolve()]
  }
];