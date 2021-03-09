
'use strict';

function assert(cond, message) {
  if (!cond) {
      message = message || "Assertion failed for unspecified reason";
      console.error(message);
      console.error(message.stack);
      throw new Exception("assert ", message);
  }
}
function _expect(cond, message) {
    if (!cond) {
        console.error(message);
    }
}
var my_assert = assert
