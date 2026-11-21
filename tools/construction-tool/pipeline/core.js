'use strict';
/* The construction tool's own geometry, for Node: the functions below are cut out of the page
   (construction-tool.template.html) verbatim, so the pipeline builds and checks every construction
   with exactly the code the browser replays it with. */
const fs = require('fs'), path = require('path');
const TPL = path.join(__dirname, '..', 'construction-tool.template.html');

const FUNCS = ['deps', 'latVec', 'polyGeom', 'curveVal', 'meet', 'evaluate', 'distClasses', 'rank', 'canon', 'cmpList', 'flexDim', 'decodeRecipe'];
const CONSTS = ['EPS', 'DTOL', 'PT', 'CURVE', 'isPt', 'CT', 'cmpc', 'sameCanon'];

function cutFunction(src, name) {
  const i = src.indexOf('function ' + name + '(');
  if (i < 0) throw new Error('core.js: function ' + name + ' not found in the template');
  let j = src.indexOf('{', i), depth = 0;
  for (let k = j; k < src.length; k++) {
    const c = src[k];
    if (c === '{') depth++;
    else if (c === '}') { depth--; if (!depth) return src.slice(i, k + 1); }
  }
  throw new Error('unbalanced ' + name);
}
function cutConst(src, name) {
  const m = src.match(new RegExp('^const ' + name + ' = .*$', 'm'));
  if (!m) throw new Error('core.js: const ' + name + ' not found');
  return m[0];
}
const src = fs.readFileSync(TPL, 'utf8');
const body = CONSTS.map(n => cutConst(src, n)).join('\n') + '\n' + FUNCS.map(n => cutFunction(src, n)).join('\n') +
  '\nreturn {' + [...CONSTS, ...FUNCS].join(',') + '};';
module.exports = new Function(body)();
