
'use strict';

function assert(cond, message) {
  if (!cond) {
      message = message || "Assertion failed for unspecified reason";
      console.error(message);
      console.error(message.stack);
      throw new Error("assert ", message);
  }
}
function _expect(cond, message) {
    if (!cond) {
        message = message || "Assertion failed for unspecified reason";
        console.error(message);
        console.error(message.stack);
    }
}
var my_assert = assert;

try {
    module.exports = {
        assert,
        _expect,
        my_assert,
        // assert3,
    };
} catch(err) {
    // not in require() mode.
}
