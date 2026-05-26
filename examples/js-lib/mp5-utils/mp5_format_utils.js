

/*
See:
  - `cModel.prototype.checkIntegrity` in the following file:
  - https://github.com/sohale/mp5-private/blob/c9f57e241c2dc4dcc584430d8483cd3cc834cbed/frontend/src/shapes/root.js

  - examples/js-lib/mp5builder.js

*/
// move to a different file
function checkMP5Object(mp5Obj) {
    const validity = IntegrityIssuesAccumulator.fine();
    if (!(mp5Obj.root && mp5Obj.root.children)) {
        validity.add(mp5Obj.root && mp5Obj.root.children, 'no .root.children', IntegrityIssuesAccumulator.HIGH_SEVERITY);
    }
    validity.react(INTEGRITY_REACTIONS.THROW_UP);
    // any reaction
    return validity;

    /*
  if (!(mp5Obj.root && mp5Obj.root.children)) throw new Error('no .root.children');
  return;
    */
}
