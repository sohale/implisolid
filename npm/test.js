'use strict';
console.log("todo: test");

var expect = require('chai').expect;

var IMPLICIT = require('./index.js');
console.log("IMPLICIT", IMPLICIT);

describe('tests', function() {
    it('initialisation (nothin)', function() {
        expect(true).to.equal(true);
        expect(!!IMPLICIT).to.equal(true);
        console.log(IMPLICIT.getLiveGeometry_from_json);
        expect(!!(IMPLICIT.getLiveGeometry_from_json)).to.equal(true);
    });
});

