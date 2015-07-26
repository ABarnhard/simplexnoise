(function() {

  function SimplexNoise(random) {
    /*
     * You can pass in a custom random number generator function if you like.
     * It is assumed to return a result between 0 & 1.
    */
    if (!!!random) { random = Math.random; }
    var i = 0;

    this.grad3 = [
      [1,1,0],[-1,1,0],[1,-1,0],[-1,-1,0],
      [1,0,1],[-1,0,1],[1,0,-1],[-1,0,-1],
      [0,1,1],[0,-1,1],[0,1,-1],[0,-1,-1]
    ];

    this.grad4 = [
      [0,1,1,1], [0,1,1,-1], [0,1,-1,1], [0,1,-1,-1],
      [0,-1,1,1], [0,-1,1,-1], [0,-1,-1,1], [0,-1,-1,-1],
      [1,0,1,1], [1,0,1,-1], [1,0,-1,1], [1,0,-1,-1],
      [-1,0,1,1], [-1,0,1,-1], [-1,0,-1,1], [-1,0,-1,-1],
      [1,1,0,1], [1,1,0,-1], [1,-1,0,1], [1,-1,0,-1],
      [-1,1,0,1], [-1,1,0,-1], [-1,-1,0,1], [-1,-1,0,-1],
      [1,1,1,0], [1,1,-1,0], [1,-1,1,0], [1,-1,-1,0],
      [-1,1,1,0], [-1,1,-1,0], [-1,-1,1,0], [-1,-1,-1,0]
    ];

    this.p = [];

    for (i = 0; i < 256; i++) {
      this.p[i] = Math.floor(random() * 256);
    }

    // To remove the need for index wrapping, double the permutation table length
    this.perm = [];
    for(i = 0; i < 512; i++) {
      this.perm[i] = this.p[i & 255];
    }

    // A lookup table to traverse the simplex around a given point in 4D.
    // Details can be found where this table is used, in the 4D noise method.
    this.simplex = [
      [0,1,2,3],[0,1,3,2],[0,0,0,0],[0,2,3,1],[0,0,0,0],[0,0,0,0],[0,0,0,0],[1,2,3,0],
      [0,2,1,3],[0,0,0,0],[0,3,1,2],[0,3,2,1],[0,0,0,0],[0,0,0,0],[0,0,0,0],[1,3,2,0],
      [0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],
      [1,2,0,3],[0,0,0,0],[1,3,0,2],[0,0,0,0],[0,0,0,0],[0,0,0,0],[2,3,0,1],[2,3,1,0],
      [1,0,2,3],[1,0,3,2],[0,0,0,0],[0,0,0,0],[0,0,0,0],[2,0,3,1],[0,0,0,0],[2,1,3,0],
      [0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],
      [2,0,1,3],[0,0,0,0],[0,0,0,0],[0,0,0,0],[3,0,1,2],[3,0,2,1],[0,0,0,0],[3,1,2,0],
      [2,1,0,3],[0,0,0,0],[0,0,0,0],[0,0,0,0],[3,1,0,2],[0,0,0,0],[3,2,0,1],[3,2,1,0]
    ];
  }

  SimplexNoise.prototype.dot = function(g, x, y) {
    return (g[0] * x) + (g[1] * y);
  };

  SimplexNoise.prototype.dot3 = function(g, x, y, z) {
    return (g[0] * x) + (g[1] * y) + (g[2] * z);
  };

  SimplexNoise.prototype.dot4 = function(g, x, y, z, w) {
    return (g[0] * x) + (g[1] * y) + (g[2] * z) + (g[3] * w);
  };

  // 2D simplex noise
  SimplexNoise.prototype.noise = function(xin, yin) {
    var n0, n1, n2, // Noise contributions from the three corners
        // Skew the input space to determine which simplex cell we're in
        F2 = 0.5 * (Math.sqrt(3) - 1),
        s  = (xin + yin) * F2, // Hairy factor for 2D
        i  = Math.floor(xin + s),
        j  = Math.floor(yin + s),
        G2 = (3 - Math.sqrt(3)) / 6,
        t  = (i + j) * G2,
        X0 = i - t, // Unskew the cell origin back to (x,y) space
        Y0 = j - t,
        x0 = xin - X0, // The x,y distances from the cell origin
        y0 = yin - Y0,
        // For the 2D case, the simplex shape is an equilateral triangle.
        // Determine which simplex we are in.
        i1, j1; // Offsets for second (middle) corner of simplex in (i,j) coords

    if (x0 > y0) {
      // lower triangle, XY order: (0,0)->(1,0)->(1,1)
      i1 = 1;
      j1 = 0;
    } else {
      // upper triangle, YX order: (0,0)->(0,1)->(1,1)
      i1 = 0;
      j1 = 1;
    }

    // A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
    // a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
    // c = (3-sqrt(3))/6
    var x1 = x0 - i1 + G2, // Offsets for middle corner in (x,y) unskewed coords
        y1 = y0 - j1 + G2,
        x2 = x0 - 1 + (2 * G2), // Offsets for last corner in (x,y) unskewed coords
        y2 = y0 - 1 + (2 * G2),
        // Work out the hashed gradient indices of the three simplex corners
        ii = i & 255,
        jj = j & 255,
        gi0 = this.perm[ii + this.perm[jj]] % 12,
        gi1 = this.perm[ii + i1 + this.perm[jj + j1]] % 12,
        gi2 = this.perm[ii + 1 + this.perm[jj + 1]] % 12,
        // Calculate the contribution from the three corners
        t0 = 0.5 - (x0 * x0) - (y0 * y0);

    if (t0 < 0) {
      n0 = 0;
    } else {
      t0 *= t0;
      n0 = t0 * t0 * this.dot(this.grad3[gi0], x0, y0);  // (x,y) of grad3 used for 2D gradient
    }

    var t1 = 0.5 - (x1 * x1) - (y1 * y1);
    if (t1 < 0) {
      n1 = 0;
    } else {
      t1 *= t1;
      n1 = t1 * t1 * this.dot(this.grad3[gi1], x1, y1);
    }

    var t2 = 0.5 - x2*x2-y2*y2;
    if (t2 < 0) {
      n2 = 0;
    } else {
      t2 *= t2;
      n2 = t2 * t2 * this.dot(this.grad3[gi2], x2, y2);
    }

    // Add contributions from each corner to get the final noise value.
    // The result is scaled to return values in the interval [-1,1].
    return 70 * (n0 + n1 + n2);
  };

  // 3D simplex noise
  SimplexNoise.prototype.noise3d = function(xin, yin, zin) {
    var n0, n1, n2, n3, // Noise contributions from the four corners
        // Skew the input space to determine which simplex cell we're in
        F3 = 1/3,
        s  = (xin + yin + zin) * F3, // Very nice and simple skew factor for 3D
        i  = Math.floor(xin + s),
        j  = Math.floor(yin + s),
        k  = Math.floor(zin + s),
        G3 = 1/6, // Very nice and simple unskew factor, too
        t  = (i + j + k) * G3,
        X0 = i - t, // Unskew the cell origin back to (x,y,z) space
        Y0 = j - t,
        Z0 = k - t,
        x0 = xin - X0, // The x,y,z distances from the cell origin
        y0 = yin - Y0,
        z0 = zin - Z0,
        // For the 3D case, the simplex shape is a slightly irregular tetrahedron.
        // Determine which simplex we are in.
        i1, j1, k1, // Offsets for second corner of simplex in (i,j,k) coords
        i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords
    if (x0 >= y0) {
      if (y0 >= z0) {
        // X Y Z order
        i1 = i2 = j2 = 1;
        j1 = k1 = k2 = 0;
      } else if (x0 >= z0) {
        // X Z Y order
        i1 = k2 = i2 = 1;
        j1 = k1 = j2 = 0;
      } else {
        // Z X Y order
        i1 = j1 = j2 = 0;
        k1 = i2 = k2 = 1;
      }
    } else {
      // x0<y0
      if (y0 < z0) {
        // Z Y X order
        i1 = j1 = i2 = 0;
        k1 = j2 = k2 = 1;
      } else if (x0 < z0) {
        // Y Z X order
        i1 = k1 = i2 = 0;
        j1 = j2 = k2 = 1;
      } else {
        // Y X Z order
        i1 = k1 = k2 = 0;
        i2 = j1 = j2 = 1;

      }
    }
    // A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
    // a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
    // a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
    // c = 1/6.
    var x1 = x0 - i1 + G3, // Offsets for second corner in (x,y,z) coords
        y1 = y0 - j1 + G3,
        z1 = z0 - k1 + G3,
        x2 = x0 - i2 + (2 * G3), // Offsets for third corner in (x,y,z) coords
        y2 = y0 - j2 + (2 * G3),
        z2 = z0 - k2 + (2 * G3),
        x3 = x0 - 1.0 + (3 * G3), // Offsets for last corner in (x,y,z) coords
        y3 = y0 - 1.0 + (3 * G3),
        z3 = z0 - 1.0 + (3 * G3),
        // Work out the hashed gradient indices of the four simplex corners
        ii = i & 255,
        jj = j & 255,
        kk = k & 255,
        gi0 = this.perm[ii + this.perm[jj + this.perm[kk]]] % 12,
        gi1 = this.perm[ii + i1 + this.perm[jj + j1 + this.perm[kk + k1]]] % 12,
        gi2 = this.perm[ii + i2 + this.perm[jj + j2 + this.perm[kk + k2]]] % 12,
        gi3 = this.perm[ii + 1 + this.perm[jj + 1 + this.perm[kk + 1]]] % 12,
        // Calculate the contribution from the four corners
        t0 = 0.6 - (x0 * x0) - (y0 * y0) - (z0 * z0);

    if (t0 < 0) {
      n0 = 0;
    } else {
      t0 *= t0;
      n0 = t0 * t0 * this.dot3(this.grad3[gi0], x0, y0, z0);
    }

    var t1 = 0.6 - (x1 * x1) - (y1 * y1) - (z1 * z1);
    if (t1 < 0) {
      n1 = 0;
    } else {
      t1 *= t1;
      n1 = t1 * t1 * this.dot3(this.grad3[gi1], x1, y1, z1);
    }

    var t2 = 0.6 - (x2 * x2) - (y2 * y2) - (z2 * z2);
    if (t2 < 0) {
      n2 = 0;
    } else {
      t2 *= t2;
      n2 = t2 * t2 * this.dot3(this.grad3[gi2], x2, y2, z2);
    }

    var t3 = 0.6 - (x3 * x3) - (y3 * y3) - (z3 * z3);
    if (t3 < 0) {
      n3 = 0;
    } else {
      t3 *= t3;
      n3 = t3 * t3 * this.dot3(this.grad3[gi3], x3, y3, z3);
    }
    // Add contributions from each corner to get the final noise value.
    // The result is scaled to stay just inside [-1,1]
    return 32 * (n0 + n1 + n2 + n3);
  };

  // 4D simplex noise
  SimplexNoise.prototype.noise4d = function( x, y, z, w ) {
    // For faster and easier lookups
    var grad4   = this.grad4,
        simplex = this.simplex,
        perm    = this.perm;

    // The skewing and unskewing factors are hairy again for the 4D case
    var F4 = (Math.sqrt(5) - 1) / 4,
        G4 = (5 - Math.sqrt(5)) / 20,
        n0, n1, n2, n3, n4, // Noise contributions from the five corners
        // Skew the (x,y,z,w) space to determine which cell of 24 simplices we're in
        s  = (x + y + z + w) * F4, // Factor for 4D skewing
        i  = Math.floor(x + s),
        j  = Math.floor(y + s),
        k  = Math.floor(z + s),
        l  = Math.floor(w + s),
        t  = (i + j + k + l) * G4, // Factor for 4D unskewing
        X0 = i - t, // Unskew the cell origin back to (x,y,z,w) space
        Y0 = j - t,
        Z0 = k - t,
        W0 = l - t,
        x0 = x - X0,  // The x,y,z,w distances from the cell origin
        y0 = y - Y0,
        z0 = z - Z0,
        w0 = w - W0;

    // For the 4D case, the simplex is a 4D shape I won't even try to describe.
    // To find out which of the 24 possible simplices we're in, we need to
    // determine the magnitude ordering of x0, y0, z0 and w0.
    // The method below is a good way of finding the ordering of x,y,z,w and
    // then find the correct traversal order for the simplex weâ€™re in.
    // First, six pair-wise comparisons are performed between each possible pair
    // of the four coordinates, and the results are used to add up binary bits
    // for an integer index.
    var c1 = (x0 > y0) ? 32 : 0,
        c2 = (x0 > z0) ? 16 : 0,
        c3 = (y0 > z0) ? 8 : 0,
        c4 = (x0 > w0) ? 4 : 0,
        c5 = (y0 > w0) ? 2 : 0,
        c6 = (z0 > w0) ? 1 : 0,
        c = c1 + c2 + c3 + c4 + c5 + c6,
        i1, j1, k1, l1, // The integer offsets for the second simplex corner
        i2, j2, k2, l2, // The integer offsets for the third simplex corner
        i3, j3, k3, l3; // The integer offsets for the fourth simplex corner
    // simplex[c] is a 4-vector with the numbers 0, 1, 2 and 3 in some order.
    // Many values of c will never occur, since e.g. x>y>z>w makes x<z, y<w and x<w
    // impossible. Only the 24 indices which have non-zero entries make any sense.
    // We use a thresholding to set the coordinates in turn from the largest magnitude.
    // The number 3 in the "simplex" array is at the position of the largest coordinate.
    i1 = simplex[c][0]>=3 ? 1 : 0;
    j1 = simplex[c][1]>=3 ? 1 : 0;
    k1 = simplex[c][2]>=3 ? 1 : 0;
    l1 = simplex[c][3]>=3 ? 1 : 0;
    // The number 2 in the "simplex" array is at the second largest coordinate.
    i2 = simplex[c][0]>=2 ? 1 : 0;
    j2 = simplex[c][1]>=2 ? 1 : 0;
    k2 = simplex[c][2]>=2 ? 1 : 0;
    l2 = simplex[c][3]>=2 ? 1 : 0;
    // The number 1 in the "simplex" array is at the second smallest coordinate.
    i3 = simplex[c][0]>=1 ? 1 : 0;
    j3 = simplex[c][1]>=1 ? 1 : 0;
    k3 = simplex[c][2]>=1 ? 1 : 0;
    l3 = simplex[c][3]>=1 ? 1 : 0;

    // The fifth corner has all coordinate offsets = 1, so no need to look that up.
    var x1  = x0 - i1 + G4, // Offsets for second corner in (x,y,z,w) coords
        y1  = y0 - j1 + G4,
        z1  = z0 - k1 + G4,
        w1  = w0 - l1 + G4,
        x2  = x0 - i2 + (2 * G4), // Offsets for third corner in (x,y,z,w) coords
        y2  = y0 - j2 + (2 * G4),
        z2  = z0 - k2 + (2 * G4),
        w2  = w0 - l2 + (2 * G4),
        x3  = x0 - i3 + (3 * G4), // Offsets for fourth corner in (x,y,z,w) coords
        y3  = y0 - j3 + (3 * G4),
        z3  = z0 - k3 + (3 * G4),
        w3  = w0 - l3 + (3 * G4),
        x4  = x0 - 1 + (4 * G4), // Offsets for last corner in (x,y,z,w) coords
        y4  = y0 - 1 + (4 * G4),
        z4  = z0 - 1 + (4 * G4),
        w4  = w0 - 1 + (4 * G4),
        // Work out the hashed gradient indices of the five simplex corners
        ii  = i & 255,
        jj  = j & 255,
        kk  = k & 255,
        ll  = l & 255,
        gi0 = perm[ii + perm[jj + perm[kk + perm[ll]]]] % 32,
        gi1 = perm[ii + i1 + perm[jj + j1 + perm[kk + k1 + perm[ll + l1]]]] % 32,
        gi2 = perm[ii + i2 + perm[jj + j2 + perm[kk + k2 + perm[ll + l2]]]] % 32,
        gi3 = perm[ii + i3 + perm[jj + j3 + perm[kk + k3 + perm[ll + l3]]]] % 32,
        gi4 = perm[ii + 1 + perm[jj + 1 + perm[kk + 1 + perm[ll + 1]]]] % 32,
        // Calculate the contribution from the five corners
        t0 = 0.6 - (x0 * x0) - (y0 * y0) - (z0 * z0) - (w0 * w0);

    if (t0 < 0) {
      n0 = 0;
    } else {
      t0 *= t0;
      n0 = t0 * t0 * this.dot4(grad4[gi0], x0, y0, z0, w0);
    }

    var t1 = 0.6 - (x1 * x1) - (y1 * y1) - (z1 * z1) - (w1 * w1);
    if (t1 < 0) {
      n1 = 0;
    } else {
      t1 *= t1;
      n1 = t1 * t1 * this.dot4(grad4[gi1], x1, y1, z1, w1);
    }

    var t2 = 0.6 - (x2 * x2) - (y2 * y2) - (z2 * z2) - (w2 * w2);
    if (t2 < 0) {
      n2 = 0;
    } else {
      t2 *= t2;
      n2 = t2 * t2 * this.dot4(grad4[gi2], x2, y2, z2, w2);
    }

    var t3 = 0.6 - (x3 * x3) - (y3 * y3) - (z3 * z3) - (w3 * w3);
    if (t3 < 0) {
      n3 = 0;
    } else {
      t3 *= t3;
      n3 = t3 * t3 * this.dot4(grad4[gi3], x3, y3, z3, w3);
    }

    var t4 = 0.6 - (x4 * x4) - (y4 * y4) - (z4 * z4) - (w4 * w4);
    if (t4 < 0) {
      n4 = 0;
    } else {
      t4 *= t4;
      n4 = t4 * t4 * this.dot4(grad4[gi4], x4, y4, z4, w4);
    }
    // Sum up and scale the result to cover the range [-1,1]
    return 27 * (n0 + n1 + n2 + n3 + n4);
  };

  // export code taken from lodash implementation

  /**
   * Checks if `value` is a global object.
   *
   * @private
   * @param {*} value The value to check.
   * @returns {null|Object} Returns `value` if it's a global object, else `null`.
   */
  function checkGlobal(value) {
    return (value && value.Object === Object) ? value : null;
  }

  /** Used to determine if values are of the language type `Object`. */
  var objectTypes = {
    'function': true,
    'object': true
  };

  /** Detect free variable `exports`. */
  var freeExports = (objectTypes[typeof exports] && exports && !exports.nodeType) ? exports : null;

  /** Detect free variable `module`. */
  var freeModule = (objectTypes[typeof module] && module && !module.nodeType) ? module : null;

  /** Detect free variable `global` from Node.js. */
  var freeGlobal = checkGlobal(freeExports && freeModule && typeof global === 'object' && global);

  /** Detect free variable `self`. */
  var freeSelf = checkGlobal(objectTypes[typeof self] && self);

  /** Detect free variable `window`. */
  var freeWindow = checkGlobal(objectTypes[typeof window] && window);

  /** Detect the popular CommonJS extension `module.exports`. */
  var moduleExports = (freeModule && freeModule.exports === freeExports) ? freeExports : null;

  /** Detect `this` as the global object. */
  var thisGlobal = checkGlobal(objectTypes[typeof this] && this);

  /**
   * Used as a reference to the global object.
   *
   * The `this` value is used if it's the global object to avoid Greasemonkey's
   * restricted `window` object, otherwise the `window` object is used.
   */
  var root = freeGlobal || ((freeWindow !== (thisGlobal && thisGlobal.window)) && freeWindow) || freeSelf || thisGlobal || Function('return this')();

  // Expose SimplexNoise on the free variable `window` or `self` when available. This
  // prevents errors in cases where SimplexNoise is loaded by a script tag in the presence
  // of an AMD loader. See http://requirejs.org/docs/errors.html#mismatch for more details.
  (freeWindow || freeSelf || {}).SimplexNoise = SimplexNoise;

  // Some AMD build optimizers like r.js check for condition patterns like the following:
  if (typeof define === 'function' && typeof define.amd === 'object' && define.amd) {
    // Define as an anonymous module so, through path mapping, it can be
    // referenced as the "underscore" module.
    define(function() {
      return SimplexNoise;
    });
  }
  // Check for `exports` after `define` in case a build optimizer adds an `exports` object.
  else if (freeExports && freeModule) {
    // Export for Node.js or RingoJS.
    if (moduleExports) {
      (freeModule.exports = SimplexNoise).SimplexNoise = SimplexNoise;
    }
    // Export for Rhino with CommonJS support.
    else {
      freeExports.SimplexNoise = SimplexNoise;
    }
  }
  else {
    // Export for a browser or Rhino.
    root.SimplexNoise = SimplexNoise;
  }
})();

