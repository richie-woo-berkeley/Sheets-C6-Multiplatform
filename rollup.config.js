import json from '@rollup/plugin-json';

export const C6_VERSION = '1.0.5';

export default {
  input: 'src/index.js',
  output: {
    file: 'dist/c6-sim.min.js',
    format: 'umd',
    name: 'C6',
    exports: 'default',
    sourcemap: true,
    banner: '// C6-Sim version 1.0.5'
  },
  plugins: [json()]
};