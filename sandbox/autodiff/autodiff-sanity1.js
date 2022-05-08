
// NodeJS code

// a fork of scripts/testing/sanity1.js

const [, , compiled_js_filename] = process.argv;
if (!compiled_js_filename) {
  throw new Error('Usage: node sanity1.js "./compiled-escripten-filename.js"');
}

async function run_main_async() {
  const wait_for_full_reload = require('../../scripts/testing/service_l1.js');
  const Service1 = await wait_for_full_reload(compiled_js_filename);
  const s1 = new Service1();
  s1.about2();

  function assert(cond, message) {
    if (!cond) {
        message = message || "Assertion failed for unspecified reason";
        console.error(message);
        console.error(message.stack);
        throw new Error("assert ", message);
    }
  }

}

run_main_async();

/*
Warning: To load an ES module, set "type": "module" in the package.json or use the .mjs extension.
*/
