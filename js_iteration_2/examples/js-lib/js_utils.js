

function __rec(obj, callback) {
  //if (Array.isArray(obj))
  callback(obj);
  for(let kk in obj) {
    callback(obj[kk]);
    __rec(obj[kk], callback);
  }
}

function freeze(obj) {
  __rec(obj, x => {Object.freeze(x);console.log('froze',x);});
  return obj;
}
