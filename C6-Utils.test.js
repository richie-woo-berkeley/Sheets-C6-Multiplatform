

import { describe, it, expect } from 'vitest';
import { merge, field, makeJSON } from 'src/C6-Utils.js';

describe('C6-Utils Utilities', () => {
  
  it('merges values with a delimiter', () => {
    const result = merge("5'", "ccata", "GGTCTC", "a", "CCAT", "cttcatccacaaatccatc", "3'", ".");
    expect(result).toBe("5'.ccata.GGTCTC.a.CCAT.cttcatccacaaatccatc.3'");
  });

  it('merges values without a delimiter error', () => {
    expect(() => merge("onlyOneArgument")).toThrow();
  });

  it('creates JSON from key-value pairs', () => {
    const input = [["name", "G00101"], ["sequence", "ATTACCGCCTTTGAGTGAGC"], ["description", "BBa_G00101 sequencing primer"]];
    const expected = '{"name":"G00101","sequence":"ATTACCGCCTTTGAGTGAGC","description":"BBa_G00101 sequencing primer"}';
    expect(makeJSON(input)).toBe(expected);
  });

  it('accesses a field from JSON', () => {
    const json = '{"name":"G00101","sequence":"ATTACCGCCTTTGAGTGAGC","description":"BBa_G00101 sequencing primer"}';
    expect(field(json, "name")).toBe("G00101");
  });

  it('throws on invalid JSON in field()', () => {
    expect(() => field("invalid-json", "name")).toThrow();
  });

});