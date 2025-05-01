import json from '@rollup/plugin-json';

export default {
  input: 'src/index.js',
  output: {
    file: 'dist/c6-sim.min.js',
    format: 'umd',
    name: 'C6',
    exports: 'default',
    sourcemap: true,
  },
  plugins: [json()]
};