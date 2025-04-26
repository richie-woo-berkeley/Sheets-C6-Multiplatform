export default {
    input: 'src/index.js',    // <-- your main code goes here
    output: {
      file: 'dist/c6-sim.min.js',
      format: 'umd',
      name: 'C6',
      exports: 'default',
      sourcemap: true
    }
  };