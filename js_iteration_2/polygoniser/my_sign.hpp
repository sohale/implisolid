#pragma once
// see implicit_vectorised_algorithms.hpp

inline REAL my_sign(REAL v, REAL ROOT_TOLERANCE) {
//     return np.sign(v) * (np.abs(v) > ROOT_TOLERANCE)
  /*
    return
      (v > 0) ?
      (+1) * (v > ROOT_TOLERANCE) :
      (v < 0) ?
      (-1) * (-v > ROOT_TOLERANCE) :
      //(v==0)
      0.0; //(0) * (v > ROOT_TOLERANCE);
  */
      /*
      REAL sgn = (v > 0) ? +1.0 : (v < 0)? -1.0 : 0.0;
      REAL vabs = (v > 0) ? v : (-v);
      // bool (v >= 0) ? (v > ROOT_TOLERANCE) : (-v > ROOT_TOLERANCE);
      REAL (vabs > ROOT_TOLERANCE) ? +1 : -1;
      return (v > 0) ? sgn
      */
    /*
    REAL r;
    if ( v > +ROOT_TOLERANCE ) {
        r = +1;
    } else if ( v < -ROOT_TOLERANCE ) {
        r = -1;
    } else {
        r = 0.0;
    }
    return r;
    */
    return
        (v > +ROOT_TOLERANCE)?
            (+1) :
        (v < -ROOT_TOLERANCE)?
            (-1)
        :
            (0.0);
}
