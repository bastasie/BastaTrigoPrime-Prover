import React, { useState, useEffect } from 'react';
import * as math from 'mathjs';
import _ from 'lodash';
import { 
  Area, AreaChart, LineChart, Line, XAxis, YAxis, CartesianGrid, 
  Tooltip, Legend, ResponsiveContainer, BarChart, Bar, ScatterChart, Scatter
} from 'recharts';

// Component for visualizing complex numbers in the complex plane
const ComplexPlane = ({ points, title }) => {
  // Transform complex points for visualization
  const transformedPoints = points.map(point => {
    const re = typeof point.re === 'function' ? point.re() : point.re;
    const im = typeof point.im === 'function' ? point.im() : point.im;
    return { x: re, y: im, name: point.name || '' };
  });

  // Find the bounds for the chart
  const maxAbs = Math.max(
    ...transformedPoints.map(p => Math.max(Math.abs(p.x), Math.abs(p.y))),
    1
  );
  const bound = Math.ceil(maxAbs * 1.2);

  return (
    <div className="bg-white p-4 rounded-lg shadow-md my-4">
      <h3 className="text-lg font-semibold mb-2">{title || "Complex Plane"}</h3>
      <ResponsiveContainer width="100%" height={300}>
        <ScatterChart
          margin={{ top: 20, right: 30, left: 20, bottom: 20 }}
        >
          <CartesianGrid strokeDasharray="3 3" />
          <XAxis 
            type="number"
            dataKey="x"
            name="Real"
            domain={[-bound, bound]}
            label={{ value: 'Re', position: 'insideBottomRight', offset: -5 }}
          />
          <YAxis 
            type="number"
            dataKey="y"
            name="Imaginary"
            domain={[-bound, bound]}
            label={{ value: 'Im', angle: -90, position: 'insideLeft' }}
          />
          <Tooltip cursor={{ strokeDasharray: '3 3' }} />
          <Scatter 
            name="Complex Numbers" 
            data={transformedPoints} 
            fill="#8884d8" 
            shape="circle"
            line={{ stroke: '#ff7300', strokeWidth: 1, strokeDasharray: '5 5' }}
          />
        </ScatterChart>
      </ResponsiveContainer>
    </div>
  );
};

// Component for visualizing polynomials
const PolynomialPlot = ({ polynomial, domain = [-10, 10], points = 100, title }) => {
  // Generate data points for the polynomial
  const data = [];
  const step = (domain[1] - domain[0]) / points;
  
  for (let x = domain[0]; x <= domain[1]; x += step) {
    try {
      const scope = { x };
      const y = math.evaluate(polynomial, scope);
      data.push({ x, y: typeof y === 'number' ? y : 0 });
    } catch (error) {
      console.error("Error evaluating polynomial:", error);
    }
  }

  return (
    <div className="bg-white p-4 rounded-lg shadow-md my-4">
      <h3 className="text-lg font-semibold mb-2">{title || `Polynomial: ${polynomial}`}</h3>
      <ResponsiveContainer width="100%" height={300}>
        <LineChart data={data} margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
          <CartesianGrid strokeDasharray="3 3" />
          <XAxis 
            dataKey="x" 
            domain={domain}
            type="number"
            label={{ value: 'x', position: 'insideBottomRight', offset: -5 }}
          />
          <YAxis 
            label={{ value: 'f(x)', angle: -90, position: 'insideLeft' }}
          />
          <Tooltip />
          <Line type="monotone" dataKey="y" stroke="#8884d8" dot={false} />
        </LineChart>
      </ResponsiveContainer>
    </div>
  );
};

// Component for visualizing matrices
const MatrixVisualizer = ({ matrix, title }) => {
  // Check if matrix is valid
  if (!matrix || !matrix.length || !matrix[0].length) {
    return <div className="text-red-500">Invalid matrix data</div>;
  }

  // Find max absolute value for color scaling
  const maxAbs = Math.max(...matrix.map(row => Math.max(...row.map(cell => Math.abs(cell)))));

  return (
    <div className="bg-white p-4 rounded-lg shadow-md my-4">
      <h3 className="text-lg font-semibold mb-2">{title || "Matrix Visualization"}</h3>
      <div className="flex justify-center">
        <div className="border-2 border-blue-500 p-2 inline-block">
          {matrix.map((row, rowIndex) => (
            <div key={rowIndex} className="flex">
              {row.map((cell, cellIndex) => (
                <div 
                  key={cellIndex} 
                  className="w-16 h-16 flex items-center justify-center border border-gray-300"
                  style={{ 
                    backgroundColor: `rgba(66, 153, 225, ${Math.min(Math.abs(cell/(maxAbs || 1)) * 0.8, 0.8)})`,
                  }}
                >
                  {typeof cell === 'number' ? cell.toFixed(2) : String(cell)}
                </div>
              ))}
            </div>
          ))}
        </div>
      </div>
    </div>
  );
};

// Component for visualizing TrigoPrime representations
const TPVisualizer = ({ trigoPrimeObj, title }) => {
  // Check if TP object is valid
  if (!trigoPrimeObj || !trigoPrimeObj.angular) {
    return <div className="text-red-500">Invalid TrigoPrime data</div>;
  }

  return (
    <div className="bg-white p-4 rounded-lg shadow-md my-4">
      <h3 className="text-lg font-semibold mb-2">{title || "TrigoPrime Representation"}</h3>
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4 p-4">
        <div className="border border-blue-300 rounded p-4">
          <h4 className="font-semibold mb-2">Angular Component</h4>
          <p className="mb-2">Value: {typeof trigoPrimeObj.angular === 'object' ? 
            `${trigoPrimeObj.angular.re.toFixed(4)} + ${trigoPrimeObj.angular.im.toFixed(4)}i` : 
            trigoPrimeObj.angular}</p>
          <p>Interpretation: {
            typeof trigoPrimeObj.angular === 'object' && trigoPrimeObj.angular.abs ?
              `Magnitude: ${trigoPrimeObj.angular.abs().toFixed(4)}, Phase: ${Math.atan2(trigoPrimeObj.angular.im, trigoPrimeObj.angular.re).toFixed(4)} radians` :
              'Cannot calculate magnitude/phase'
          }</p>
        </div>
        <div className="border border-green-300 rounded p-4">
          <h4 className="font-semibold mb-2">Prime Component</h4>
          <p className="mb-2">Value: {typeof trigoPrimeObj.prime === 'object' ? 
            JSON.stringify(trigoPrimeObj.prime) : 
            String(trigoPrimeObj.prime)}</p>
          <p>Encoding: {typeof trigoPrimeObj.primeFactors === 'object' ? 
            Object.entries(trigoPrimeObj.primeFactors).map(([p, exp]) => `${p}^${exp}`).join(' · ') : 
            'Prime factorization not available'}</p>
        </div>
      </div>
      <div className="mt-4 p-4 border border-purple-300 rounded">
        <h4 className="font-semibold mb-2">Original Value Decoded</h4>
        <p>Value: {typeof trigoPrimeObj.original === 'object' ? 
          JSON.stringify(trigoPrimeObj.original) : 
          String(trigoPrimeObj.original)}</p>
      </div>
    </div>
  );
};

// Component for visualizing Fourier Series
const FourierSeriesVisualizer = ({ originalFn, terms = 5, domain = [-Math.PI, Math.PI], points = 200, title }) => {
  const [seriesData, setSeriesData] = useState([]);
  const [coefficients, setCoefficients] = useState({ a0: 0, a: [], b: [] });
  
  useEffect(() => {
    const data = [];
    const step = (domain[1] - domain[0]) / points;
    
    // Generate coefficients
    const a = [];
    const b = [];
    
    // For simplicity, we'll approximate the integrals
    const n = terms;
    const dx = (domain[1] - domain[0]) / 1000;
    
    // Calculate a0
    let a0 = 0;
    for (let x = domain[0]; x <= domain[1]; x += dx) {
      try {
        const scope = { x };
        const fx = math.evaluate(originalFn, scope);
        a0 += fx * dx;
      } catch (error) {
        console.error("Error calculating a0:", error);
      }
    }
    a0 = a0 / (domain[1] - domain[0]);
    
    // Calculate a_n and b_n for n=1 to terms
    for (let k = 1; k <= n; k++) {
      let ak = 0;
      let bk = 0;
      
      for (let x = domain[0]; x <= domain[1]; x += dx) {
        try {
          const scope = { x };
          const fx = math.evaluate(originalFn, scope);
          const cosVal = Math.cos((2 * Math.PI * k * (x - domain[0])) / (domain[1] - domain[0]));
          const sinVal = Math.sin((2 * Math.PI * k * (x - domain[0])) / (domain[1] - domain[0]));
          
          ak += fx * cosVal * dx;
          bk += fx * sinVal * dx;
        } catch (error) {
          console.error(`Error calculating a${k} or b${k}:`, error);
        }
      }
      
      ak = (2 * ak) / (domain[1] - domain[0]);
      bk = (2 * bk) / (domain[1] - domain[0]);
      
      a.push(ak);
      b.push(bk);
    }
    
    setCoefficients({ a0, a, b });
    
    // Generate data for plotting
    for (let x = domain[0]; x <= domain[1]; x += step) {
      try {
        const scope = { x };
        const originalY = math.evaluate(originalFn, scope);
        
        // Calculate Fourier series approximation
        let fourierY = a0 / 2;
        for (let k = 0; k < n; k++) {
          const cosVal = Math.cos((2 * Math.PI * (k + 1) * (x - domain[0])) / (domain[1] - domain[0]));
          const sinVal = Math.sin((2 * Math.PI * (k + 1) * (x - domain[0])) / (domain[1] - domain[0]));
          
          fourierY += a[k] * cosVal + b[k] * sinVal;
        }
        
        data.push({ x, original: originalY, fourier: fourierY });
      } catch (error) {
        console.error("Error calculating series values:", error);
      }
    }
    
    setSeriesData(data);
  }, [originalFn, terms, domain, points]);
  
  return (
    <div className="bg-white p-4 rounded-lg shadow-md my-4">
      <h3 className="text-lg font-semibold mb-2">{title || `Fourier Series Approximation: ${terms} terms`}</h3>
      <div className="mb-4">
        <h4 className="font-medium mb-2">Fourier Coefficients:</h4>
        <p className="mb-1">a₀ = {coefficients.a0.toFixed(4)}</p>
        <div className="grid grid-cols-2 gap-2">
          <div>
            {coefficients.a.map((val, idx) => (
              <p key={`a${idx+1}`} className="mb-1">a₍{idx+1}₎ = {val.toFixed(4)}</p>
            ))}
          </div>
          <div>
            {coefficients.b.map((val, idx) => (
              <p key={`b${idx+1}`} className="mb-1">b₍{idx+1}₎ = {val.toFixed(4)}</p>
            ))}
          </div>
        </div>
      </div>
      <ResponsiveContainer width="100%" height={300}>
        <LineChart data={seriesData} margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
          <CartesianGrid strokeDasharray="3 3" />
          <XAxis 
            dataKey="x" 
            domain={domain}
            type="number"
            label={{ value: 'x', position: 'insideBottomRight', offset: -5 }}
          />
          <YAxis />
          <Tooltip />
          <Legend />
          <Line type="monotone" dataKey="original" stroke="#8884d8" name="Original Function" dot={false} />
          <Line type="monotone" dataKey="fourier" stroke="#82ca9d" name="Fourier Approximation" dot={false} />
        </LineChart>
      </ResponsiveContainer>
    </div>
  );
};

// Component for visualizing wave functions in quantum mechanics
const WaveFunctionVisualizer = ({ psiRe, psiIm, domain = [-10, 10], points = 200, title }) => {
  // Generate data points for the wave function
  const generateData = () => {
    const data = [];
    const step = (domain[1] - domain[0]) / points;
    
    for (let x = domain[0]; x <= domain[1]; x += step) {
      try {
        const scope = { x };
        const real = math.evaluate(psiRe, scope);
        const imag = math.evaluate(psiIm, scope);
        const probability = real**2 + imag**2;  // |ψ|²
        
        data.push({ x, real, imag, probability });
      } catch (error) {
        console.error("Error evaluating wave function:", error);
      }
    }
    
    return data;
  };
  
  const waveData = generateData();
  
  return (
    <div className="bg-white p-4 rounded-lg shadow-md my-4">
      <h3 className="text-lg font-semibold mb-2">{title || "Quantum Wave Function"}</h3>
      <div className="mb-4">
        <p className="font-medium">Wave Function: ψ(x) = {psiRe} + i·({psiIm})</p>
        <p className="text-sm text-gray-600 mt-1">The probability density |ψ|² represents the likelihood of finding a particle at position x</p>
      </div>
      <ResponsiveContainer width="100%" height={300}>
        <LineChart data={waveData} margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
          <CartesianGrid strokeDasharray="3 3" />
          <XAxis 
            dataKey="x" 
            domain={domain}
            type="number"
            label={{ value: 'x', position: 'insideBottomRight', offset: -5 }}
          />
          <YAxis />
          <Tooltip />
          <Legend />
          <Line type="monotone" dataKey="real" stroke="#8884d8" name="Re(ψ)" dot={false} />
          <Line type="monotone" dataKey="imag" stroke="#82ca9d" name="Im(ψ)" dot={false} />
          <Line type="monotone" dataKey="probability" stroke="#ff7300" name="|ψ|²" dot={false} />
        </LineChart>
      </ResponsiveContainer>
    </div>
  );
};

// Main TrigoPrime class for mathematical operations
class TrigoPrime {
  // Prime numbers cache for optimization
  static primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];
  
  // Get the nth prime number
  static getNthPrime(n) {
    if (n <= this.primes.length) {
      return this.primes[n-1];
    }
    
    // If we need more primes, generate them (simple method for demonstration)
    let candidate = this.primes[this.primes.length - 1] + 2;
    while (this.primes.length < n) {
      let isPrime = true;
      for (let i = 0; i < this.primes.length && this.primes[i] * this.primes[i] <= candidate; i++) {
        if (candidate % this.primes[i] === 0) {
          isPrime = false;
          break;
        }
      }
      if (isPrime) {
        this.primes.push(candidate);
      }
      candidate += 2;
    }
    
    return this.primes[n-1];
  }
  
  // Prime factorization of a number
  static primeFactorize(num) {
    if (num <= 1) return { 1: 1 };
    
    const absNum = Math.abs(num);
    const factors = {};
    
    // Special case for negative numbers in our encoding
    if (num < 0) {
      factors[2] = 1;
      num = absNum;
    } else {
      num = absNum;
    }
    
    // Trial division
    for (let i = 0; i < this.primes.length && this.primes[i] * this.primes[i] <= num; i++) {
      const p = this.primes[i];
      while (num % p === 0) {
        factors[p] = (factors[p] || 0) + 1;
        num /= p;
      }
    }
    
    // If num > 1, it's a prime number itself
    if (num > 1) {
      factors[num] = (factors[num] || 0) + 1;
    }
    
    return factors;
  }
  
  // Basic conversion functions
  static toTrigoRepresentation(x) {
    if (typeof x === 'number') {
      // For integers
      if (Number.isInteger(x)) {
        // Angular component: e^(iπx)
        const angular = math.complex(Math.cos(Math.PI * x), Math.sin(Math.PI * x));
        
        // Prime encoding for integers based on the TrigoPrime paper
        let prime;
        if (x === 0) {
          prime = 1;  // TP(0) = {1, 1}
        } else if (x > 0) {
          prime = this.getNthPrime(x);  // TP(n) = {e^(iπn), prime(n)}
        } else {
          prime = 2 * this.getNthPrime(-x);  // TP(-n) = {e^(-iπn), 2·prime(n)}
        }
        
        // Get prime factorization for display purposes
        const primeFactors = this.primeFactorize(x);
        
        return { 
          angular, 
          prime, 
          primeFactors,
          original: x,
          type: 'integer'
        };
      }
      // For rational numbers
      else if (isFinite(x)) {
        const [num, den] = this.toFraction(x);
        
        // Angular component: e^(iπ·num/den)
        const angular = math.complex(
          Math.cos(Math.PI * (num/den)), 
          Math.sin(Math.PI * (num/den))
        );
        
        // Prime encodings for numerator and denominator
        const primeNum = num === 0 ? 1 : (num > 0 ? this.getNthPrime(num) : 2 * this.getNthPrime(-num));
        const primeDen = den === 0 ? 1 : (den > 0 ? this.getNthPrime(den) : 2 * this.getNthPrime(-den));
        
        return { 
          angular, 
          prime: primeNum / primeDen, 
          primeFactors: { numerator: this.primeFactorize(num), denominator: this.primeFactorize(den) },
          original: x,
          type: 'rational',
          rational: { numerator: num, denominator: den }
        };
      }
    }
    // For complex numbers
    else if (x && typeof x === 'object' && x.re !== undefined && x.im !== undefined) {
      const re = typeof x.re === 'function' ? x.re() : x.re;
      const im = typeof x.im === 'function' ? x.im() : x.im;
      
      const magnitude = Math.sqrt(re**2 + im**2);
      const angle = Math.atan2(im, re);
      
      // Angular component: e^(iπ|z|e^(i·arg(z)))
      const innerComplex = math.complex(
        magnitude * Math.cos(angle),
        magnitude * Math.sin(angle)
      );
      
      const angular = math.complex(
        Math.cos(Math.PI * innerComplex.re),
        Math.sin(Math.PI * innerComplex.im)
      );
      
      // Prime encoding for complex numbers (simplified for demonstration)
      // In full implementation, would use a more sophisticated encoding
      const primeRe = Number.isInteger(re) ? 
        (re === 0 ? 1 : (re > 0 ? this.getNthPrime(re) : 2 * this.getNthPrime(-re))) : 
        `E(${re})`;
      
      const primeIm = Number.isInteger(im) ? 
        (im === 0 ? 1 : (im > 0 ? this.getNthPrime(im) : 2 * this.getNthPrime(-im))) : 
        `E(${im})`;
      
      return { 
        angular, 
        prime: { real: primeRe, imaginary: primeIm }, 
        primeFactors: { real: this.primeFactorize(re), imaginary: this.primeFactorize(im) },
        original: { re, im },
        type: 'complex'
      };
    }
    
    // Default fallback
    return { 
      angular: null, 
      prime: null, 
      primeFactors: null,
      original: x,
      type: 'unknown'
    };
  }
  
  // Decode TrigoPrime representation back to original value
  static fromTrigoRepresentation(tpObj) {
    if (!tpObj) return null;
    
    // If we have the original value, return it directly
    if (tpObj.original !== undefined) {
      return tpObj.original;
    }
    
    // Try to decode based on the type
    if (tpObj.type === 'integer') {
      // For integers, determine the value from the prime encoding
      // In a real implementation, this would have a more robust algorithm
      if (tpObj.prime === 1) return 0;
      
      // Check if it's a positive or negative number
      const primeIdx = this.primes.indexOf(tpObj.prime);
      if (primeIdx >= 0) {
        return primeIdx + 1;  // Positive integer
      }
      
      // Check if it's a negative integer (prime = 2 * prime(n))
      if (tpObj.prime % 2 === 0) {
        const halfPrime = tpObj.prime / 2;
        const primeIdx = this.primes.indexOf(halfPrime);
        if (primeIdx >= 0) {
          return -(primeIdx + 1);  // Negative integer
        }
      }
    } 
    else if (tpObj.type === 'rational') {
      if (tpObj.rational) {
        return tpObj.rational.numerator / tpObj.rational.denominator;
      }
    }
    else if (tpObj.type === 'complex') {
      if (tpObj.original) {
        return math.complex(tpObj.original.re, tpObj.original.im);
      }
    }
    
    // If we couldn't decode, return the original representation
    return tpObj;
  }
  
  // Convert to fraction representation
  static toFraction(decimal) {
    const tolerance = 1.0E-10;
    let h1 = 1;
    let h2 = 0;
    let k1 = 0;
    let k2 = 1;
    let b = decimal;
    
    do {
      let a = Math.floor(b);
      let aux = h1;
      h1 = a * h1 + h2;
      h2 = aux;
      aux = k1;
      k1 = a * k1 + k2;
      k2 = aux;
      b = 1 / (b - a);
    } while (Math.abs(decimal - h1 / k1) > decimal * tolerance && isFinite(b));
    
    return [h1, k1];
  }
  
  // Addition in TrigoPrime
  static add(a, b) {
    const aTP = this.toTrigoRepresentation(a);
    const bTP = this.toTrigoRepresentation(b);
    
    // Direct addition for the original values
    const result = typeof aTP.original === 'number' && typeof bTP.original === 'number' ?
      aTP.original + bTP.original :
      null;
    
    // TP representation of the result
    const resultTP = result !== null ? 
      this.toTrigoRepresentation(result) : 
      {
        // Calculate result components based on TrigoPrime rules
        angular: math.complex(
          Math.cos(Math.PI * (typeof aTP.original === 'number' && typeof bTP.original === 'number' ? 
            aTP.original + bTP.original : 0)),
          Math.sin(Math.PI * (typeof aTP.original === 'number' && typeof bTP.original === 'number' ? 
            aTP.original + bTP.original : 0))
        ),
        prime: `E(${aTP.original} + ${bTP.original})`,
        original: { operation: 'add', operands: [aTP.original, bTP.original] },
        type: 'expression'
      };
    
    // Decoded value
    const decoded = this.fromTrigoRepresentation(resultTP);
    
    return {
      result: result !== null ? result : `${aTP.original} + ${bTP.original}`,
      tpResult: resultTP,
      decoded,
      operation: 'addition',
      operands: [aTP, bTP]
    };
  }
  
  // Multiplication in TrigoPrime
  static multiply(a, b) {
    const aTP = this.toTrigoRepresentation(a);
    const bTP = this.toTrigoRepresentation(b);
    
    // Direct multiplication for the original values
    const result = typeof aTP.original === 'number' && typeof bTP.original === 'number' ?
      aTP.original * bTP.original :
      null;
    
    // TP representation of the result
    // For multiplication: TP(a·b) = {e^(iπa)·e^(iπb), E(a)·E(b)}
    const resultTP = result !== null ? 
      this.toTrigoRepresentation(result) : 
      {
        // Calculate result components based on TrigoPrime rules
        angular: aTP.angular && bTP.angular ? 
          math.multiply(aTP.angular, bTP.angular) : 
          null,
        prime: typeof aTP.prime === 'number' && typeof bTP.prime === 'number' ? 
          aTP.prime * bTP.prime : 
          `${aTP.prime} · ${bTP.prime}`,
        original: { operation: 'multiply', operands: [aTP.original, bTP.original] },
        type: 'expression'
      };
    
    // Decoded value
    const decoded = this.fromTrigoRepresentation(resultTP);
    
    return {
      result: result !== null ? result : `${aTP.original} · ${bTP.original}`,
      tpResult: resultTP,
      decoded,
      operation: 'multiplication',
      operands: [aTP, bTP]
    };
  }
  
  // Matrix operations
  static matrixAdd(A, B) {
    if (!A || !B || A.length !== B.length || A[0].length !== B.length) {
      throw new Error("Matrices must have the same dimensions for addition");
    }
    
    const result = [];
    const tpMatrixA = [];
    const tpMatrixB = [];
    const tpResult = [];
    
    // Convert input matrices to TrigoPrime representations
    for (let i = 0; i < A.length; i++) {
      tpMatrixA[i] = [];
      tpMatrixB[i] = [];
      for (let j = 0; j < A[0].length; j++) {
        tpMatrixA[i][j] = this.toTrigoRepresentation(A[i][j]);
        tpMatrixB[i][j] = this.toTrigoRepresentation(B[i][j]);
      }
    }
    
    // Perform matrix addition
    for (let i = 0; i < A.length; i++) {
      result[i] = [];
      tpResult[i] = [];
      for (let j = 0; j < A[0].length; j++) {
        result[i][j] = A[i][j] + B[i][j];
        tpResult[i][j] = this.toTrigoRepresentation(result[i][j]);
      }
    }
    
    // Provide a decoded version
    const decoded = result.map(row => row.slice());
    
    return {
      result,
      tpMatrixA,
      tpMatrixB,
      tpResult,
      decoded,
      operation: 'matrix-addition',
      operands: [A, B]
    };
  }
  
  static matrixMultiply(A, B) {
    if (!A || !B || A[0].length !== B.length) {
      throw new Error("Matrix dimensions are not compatible for multiplication");
    }
    
    const result = [];
    const tpMatrixA = [];
    const tpMatrixB = [];
    const tpResult = [];
    
    // Convert input matrices to TrigoPrime representations
    for (let i = 0; i < A.length; i++) {
      tpMatrixA[i] = [];
      for (let j = 0; j < A[0].length; j++) {
        tpMatrixA[i][j] = this.toTrigoRepresentation(A[i][j]);
      }
    }
    
    for (let i = 0; i < B.length; i++) {
      tpMatrixB[i] = [];
      for (let j = 0; j < B[0].length; j++) {
        tpMatrixB[i][j] = this.toTrigoRepresentation(B[i][j]);
      }
    }
    
    // Perform matrix multiplication
    for (let i = 0; i < A.length; i++) {
      result[i] = [];
      tpResult[i] = [];
      for (let j = 0; j < B[0].length; j++) {
        result[i][j] = 0;
        for (let k = 0; k < A[0].length; k++) {
          result[i][j] += A[i][k] * B[k][j];
        }
        tpResult[i][j] = this.toTrigoRepresentation(result[i][j]);
      }
    }
    
    // Provide a decoded version
    const decoded = result.map(row => row.slice());
    
    return {
      result,
      tpMatrixA,
      tpMatrixB,
      tpResult,
      decoded,
      operation: 'matrix-multiplication',
      operands: [A, B]
    };
  }
  
  // Calculus operations
  static derivative(expression, variable = 'x') {
    try {
      const derivative = math.derivative(expression, variable).toString();
      
      // Generate sample points to visualize the derivative
      const domain = [-5, 5];
      const points = 100;
      const step = (domain[1] - domain[0]) / points;
      
      const data = [];
      for (let x = domain[0]; x <= domain[1]; x += step) {
        try {
          const scope = { [variable]: x };
          const original = math.evaluate(expression, scope);
          const derived = math.evaluate(derivative, scope);
          data.push({ x, original, derived });
        } catch (error) {
          console.error("Error evaluating functions:", error);
        }
      }
      
      // TrigoPrime representation
      const tpExpression = {
        angular: `e^(iπ${expression})`,
        prime: `E(${expression})`,
        original: expression,
        type: 'expression'
      };
      
      const tpDerivative = {
        angular: `e^(iπ${derivative})`,
        prime: `E(${derivative})`,
        original: derivative,
        type: 'expression'
      };
      
      return {
        result: derivative,
        data,
        tpExpression,
        tpDerivative,
        decoded: derivative,
        operation: 'derivative',
        operand: expression
      };
    } catch (error) {
      console.error("Error calculating derivative:", error);
      return {
        error: error.message,
        operation: 'derivative',
        operand: expression
      };
    }
  }
  
  // Solving differential equations (simplified)
  static solveDifferentialEq(equation, initialConditions = {}) {
    // For demonstration purposes, we'll solve simple first-order ODEs
    // In a full implementation, this would be much more complex
    
    // Generate numerical solution for visualization
    const numerical = this.generateNumericalSolution(equation, initialConditions);
    
    // For this demo, we'll use a placeholder analytical solution
    // In a real implementation, this would use proper differential equation solvers
    const analyticalSolution = equation === 'y' ? 
      `y = C·e^x` : 
      `y = [Analytical solution placeholder]`;
    
    // TrigoPrime representation
    const tpEquation = {
      angular: `e^(iπ(${equation}))`,
      prime: `E(${equation})`,
      original: equation,
      type: 'expression'
    };
    
    const tpSolution = {
      angular: `e^(iπ(${analyticalSolution}))`,
      prime: `E(${analyticalSolution})`,
      original: analyticalSolution,
      type: 'expression'
    };
    
    return {
      result: analyticalSolution,
      numerical,
      tpEquation,
      tpSolution,
      decoded: analyticalSolution,
      operation: 'differential-equation',
      operand: equation
    };
  }
  
  // Generate a simple numerical solution for visualization
  static generateNumericalSolution(equation, initialConditions) {
    // This is a placeholder that would be replaced with actual numerical methods
    const data = [];
    const x0 = initialConditions.x || 0;
    const y0 = initialConditions.y || 1;
    
    let x = x0;
    let y = y0;
    const h = 0.1; // Step size
    
    // Simple Euler method for demonstration
    for (let i = 0; i < 100; i++) {
      data.push({ x, y });
      
      // Evaluate the right-hand side of dy/dx = f(x,y)
      let slope;
      try {
        // Attempt to evaluate the equation
        const scope = { x, y };
        slope = math.evaluate(equation, scope);
      } catch (error) {
        // Fallback to exponential solution
        slope = y;
      }
      
      y = y + h * slope;
      x = x + h;
    }
    
    return data;
  }
  
  // Fourier series calculations
  static computeFourierSeries(expression, terms = 5, domain = [-Math.PI, Math.PI]) {
    try {
      // For simplicity, we'll approximate the integrals
      const dx = (domain[1] - domain[0]) / 1000;
      
      // Calculate a0
      let a0 = 0;
      for (let x = domain[0]; x <= domain[1]; x += dx) {
        try {
          const scope = { x };
          const fx = math.evaluate(expression, scope);
          a0 += fx * dx;
        } catch (error) {
          console.error("Error calculating a0:", error);
        }
      }
      a0 = a0 / (domain[1] - domain[0]);
      
      // Calculate a_n and b_n for n=1 to terms
      const a = [];
      const b = [];
      
      for (let k = 1; k <= terms; k++) {
        let ak = 0;
        let bk = 0;
        
        for (let x = domain[0]; x <= domain[1]; x += dx) {
          try {
            const scope = { x };
            const fx = math.evaluate(expression, scope);
            const cosVal = Math.cos((2 * Math.PI * k * (x - domain[0])) / (domain[1] - domain[0]));
            const sinVal = Math.sin((2 * Math.PI * k * (x - domain[0])) / (domain[1] - domain[0]));
            
            ak += fx * cosVal * dx;
            bk += fx * sinVal * dx;
          } catch (error) {
            console.error(`Error calculating a${k} or b${k}:`, error);
          }
        }
        
        ak = (2 * ak) / (domain[1] - domain[0]);
        bk = (2 * bk) / (domain[1] - domain[0]);
        
        a.push(ak);
        b.push(bk);
      }
      
      // Generate the Fourier series expression
      let fourierExpression = `${(a0/2).toFixed(4)}`;
      for (let k = 0; k < terms; k++) {
        if (Math.abs(a[k]) > 1e-10) {
          fourierExpression += ` + ${a[k].toFixed(4)} * cos(${k+1}x)`;
        }
        if (Math.abs(b[k]) > 1e-10) {
          fourierExpression += ` + ${b[k].toFixed(4)} * sin(${k+1}x)`;
        }
      }
      
      // TrigoPrime representations
      const tpOriginalFunction = {
        angular: `e^(iπ(${expression}))`,
        prime: `E(${expression})`,
        original: expression,
        type: 'expression'
      };
      
      const tpFourierSeries = {
        angular: `e^(iπ(${fourierExpression}))`,
        prime: `E(${fourierExpression})`,
        original: fourierExpression,
        type: 'expression'
      };
      
      return {
        originalFunction: expression,
        fourierSeries: fourierExpression,
        coefficients: { a0, a, b },
        tpOriginalFunction,
        tpFourierSeries,
        decoded: fourierExpression,
        operation: 'fourier-series',
        terms
      };
    } catch (error) {
      console.error("Error computing Fourier series:", error);
      return {
        error: error.message,
        operation: 'fourier-series',
        operand: expression
      };
    }
  }
  
  // Quantum wave function analysis
  static analyzeWaveFunction(psiRe, psiIm) {
    try {
      // Full wave function expression
      const waveFunction = `(${psiRe}) + i*(${psiIm})`;
      
      // TrigoPrime representation
      const tpWaveFunction = {
        angular: `e^(iπ|ψ|e^(i*arg(ψ)))`,
        prime: `E(${waveFunction})`,
        original: waveFunction,
        type: 'wave-function'
      };
      
      // Generate characteristic data about the wave function
      const domain = [-10, 10];
      const points = 100;
      const step = (domain[1] - domain[0]) / points;
      
      let maxProbability = 0;
      let maxProbabilityPosition = 0;
      let totalProbability = 0;
      
      for (let x = domain[0]; x <= domain[1]; x += step) {
        try {
          const scope = { x };
          const real = math.evaluate(psiRe, scope);
          const imag = math.evaluate(psiIm, scope);
          const probability = real**2 + imag**2;
          
          totalProbability += probability * step;
          
          if (probability > maxProbability) {
            maxProbability = probability;
            maxProbabilityPosition = x;
          }
        } catch (error) {
          console.error("Error evaluating wave function:", error);
        }
      }
      
      // Normalization factor (for a properly normalized wave function, this would be 1)
      const normalizationFactor = 1 / Math.sqrt(totalProbability);
      
      return {
        waveFunction,
        psiRe,
        psiIm,
        maxProbability,
        maxProbabilityPosition,
        totalProbability,
        normalizationFactor,
        tpWaveFunction,
        decoded: waveFunction,
        operation: 'wave-function-analysis'
      };
    } catch (error) {
      console.error("Error analyzing wave function:", error);
      return {
        error: error.message,
        operation: 'wave-function-analysis',
        operands: [psiRe, psiIm]
      };
    }
  }
}

// Main application component
function TrigoprimeApp() {
  const [activeTab, setActiveTab] = useState('arithmetic');
  const [inputA, setInputA] = useState('2');
  const [inputB, setInputB] = useState('3');
  const [matrixA, setMatrixA] = useState([[1, 2], [3, 4]]);
  const [matrixB, setMatrixB] = useState([[5, 6], [7, 8]]);
  const [polynomialExpr, setPolynomialExpr] = useState('2*x^2 + 3*x + 1');
  const [calcExpr, setCalcExpr] = useState('x^2');
  const [diffEqExpr, setDiffEqExpr] = useState('y');
  const [fourierExpr, setFourierExpr] = useState('x');
  const [waveReExpr, setWaveReExpr] = useState('cos(x)');
  const [waveImExpr, setWaveImExpr] = useState('sin(x)');
  const [fourierTerms, setFourierTerms] = useState(5);
  const [result, setResult] = useState(null);
  
  // Handle form submissions
  const handleArithmeticSubmit = (operation) => {
    try {
      const a = parseFloat(inputA);
      const b = parseFloat(inputB);
      
      if (isNaN(a) || isNaN(b)) {
        throw new Error("Invalid input: Please enter valid numbers");
      }
      
      let res;
      if (operation === 'add') {
        res = TrigoPrime.add(a, b);
      } else if (operation === 'multiply') {
        res = TrigoPrime.multiply(a, b);
      }
      
      setResult(res);
    } catch (error) {
      console.error("Error:", error);
      setResult({ error: error.message });
    }
  };
  
  const handleMatrixSubmit = (operation) => {
    try {
      let res;
      if (operation === 'add') {
        res = TrigoPrime.matrixAdd(matrixA, matrixB);
      } else if (operation === 'multiply') {
        res = TrigoPrime.matrixMultiply(matrixA, matrixB);
      }
      
      setResult(res);
    } catch (error) {
      console.error("Error:", error);
      setResult({ error: error.message });
    }
  };
  
  const handlePolynomialSubmit = () => {
    try {
      // For polynomials, we'll just visualize them
      setResult({
        operation: 'polynomial',
        expression: polynomialExpr,
        tpExpression: {
          angular: `e^(iπ(${polynomialExpr}))`,
          prime: `E(${polynomialExpr})`,
          original: polynomialExpr,
          type: 'expression'
        },
        decoded: polynomialExpr
      });
    } catch (error) {
      console.error("Error:", error);
      setResult({ error: error.message });
    }
  };
  
  const handleCalculusSubmit = () => {
    try {
      const res = TrigoPrime.derivative(calcExpr);
      setResult(res);
    } catch (error) {
      console.error("Error:", error);
      setResult({ error: error.message });
    }
  };
  
  const handleDiffEqSubmit = () => {
    try {
      const res = TrigoPrime.solveDifferentialEq(diffEqExpr, { x: 0, y: 1 });
      setResult(res);
    } catch (error) {
      console.error("Error:", error);
      setResult({ error: error.message });
    }
  };
  
  const handleFourierSubmit = () => {
    try {
      const res = TrigoPrime.computeFourierSeries(fourierExpr, fourierTerms);
      setResult(res);
    } catch (error) {
      console.error("Error:", error);
      setResult({ error: error.message });
    }
  };
  
  const handleWaveFunctionSubmit = () => {
    try {
      const res = TrigoPrime.analyzeWaveFunction(waveReExpr, waveImExpr);
      setResult(res);
    } catch (error) {
      console.error("Error:", error);
      setResult({ error: error.message });
    }
  };
  
  // Handle matrix input changes
  const handleMatrixCellChange = (matrix, setter, i, j, value) => {
    const newMatrix = [...matrix];
    newMatrix[i][j] = parseFloat(value) || 0;
    setter(newMatrix);
  };
  
  // Render result visualization based on the operation
  const renderResult = () => {
    if (!result) return null;
    
    if (result.error) {
      return <div className="text-red-500 p-4">{result.error}</div>;
    }
    
    switch (result.operation) {
      case 'addition':
      case 'multiplication':
        return (
          <div className="bg-white p-6 rounded-lg shadow-md">
            <h3 className="text-xl font-semibold mb-4">Result: {result.result}</h3>
            
            <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
              <div className="border p-4 rounded bg-gray-50">
                <h4 className="font-medium mb-2">Input A:</h4>
                <p>Value: {result.operands[0].original}</p>
                <p>TP Representation: {`{${typeof result.operands[0].angular === 'object' ? 
                  `${result.operands[0].angular.re.toFixed(4)} + ${result.operands[0].angular.im.toFixed(4)}i` : 
                  result.operands[0].angular}, ${result.operands[0].prime}}`}</p>
              </div>
              
              <div className="border p-4 rounded bg-gray-50">
                <h4 className="font-medium mb-2">Input B:</h4>
                <p>Value: {result.operands[1].original}</p>
                <p>TP Representation: {`{${typeof result.operands[1].angular === 'object' ? 
                  `${result.operands[1].angular.re.toFixed(4)} + ${result.operands[1].angular.im.toFixed(4)}i` : 
                  result.operands[1].angular}, ${result.operands[1].prime}}`}</p>
              </div>
            </div>
            
            <div className="border p-4 rounded bg-blue-50 mb-4">
              <h4 className="font-medium mb-2">TrigoPrime Representation of Result:</h4>
              <p>Angular Component: {typeof result.tpResult.angular === 'object' ? 
                `${result.tpResult.angular.re.toFixed(4)} + ${result.tpResult.angular.im.toFixed(4)}i` : 
                result.tpResult.angular}</p>
              <p>Prime Component: {typeof result.tpResult.prime === 'object' ? 
                JSON.stringify(result.tpResult.prime) : 
                result.tpResult.prime}</p>
            </div>
            
            <div className="border p-4 rounded bg-green-50">
              <h4 className="font-medium mb-2">Decoded Result:</h4>
              <p>{typeof result.decoded === 'object' ? 
                JSON.stringify(result.decoded) : 
                result.decoded}</p>
            </div>
          </div>
        );
      
      case 'matrix-addition':
      case 'matrix-multiplication':
        return (
          <div className="bg-white p-6 rounded-lg shadow-md">
            <h3 className="text-xl font-semibold mb-4">Matrix {result.operation === 'matrix-addition' ? 'Addition' : 'Multiplication'} Result:</h3>
            
            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
              <MatrixVisualizer matrix={result.operands[0]} title="Matrix A" />
              <MatrixVisualizer matrix={result.operands[1]} title="Matrix B" />
              <div className="md:col-span-2">
                <MatrixVisualizer matrix={result.result} title="Result Matrix" />
              </div>
            </div>
            
            <div className="mt-6 border p-4 rounded bg-blue-50">
              <h4 className="font-medium mb-2">TrigoPrime Matrix Representation:</h4>
              <p className="mb-2">Each cell in the matrix is encoded using the TrigoPrime dual representation.</p>
              <p>For example, cell (0,0): {typeof result.tpResult[0][0].angular === 'object' ? 
                `{${result.tpResult[0][0].angular.re.toFixed(4)} + ${result.tpResult[0][0].angular.im.toFixed(4)}i, ${result.tpResult[0][0].prime}}` : 
                `{${result.tpResult[0][0].angular}, ${result.tpResult[0][0].prime}}`}</p>
            </div>
            
            <div className="mt-4 border p-4 rounded bg-green-50">
              <h4 className="font-medium mb-2">Decoded Result Matrix:</h4>
              <pre className="overflow-x-auto">{JSON.stringify(result.decoded, null, 2)}</pre>
            </div>
          </div>
        );
      
      case 'polynomial':
        return (
          <div className="bg-white p-6 rounded-lg shadow-md">
            <h3 className="text-xl font-semibold mb-4">Polynomial Visualization:</h3>
            <PolynomialPlot polynomial={result.expression} domain={[-5, 5]} title={`P(x) = ${result.expression}`} />
            
            <div className="mt-6 border p-4 rounded bg-blue-50">
              <h4 className="font-medium mb-2">TrigoPrime Polynomial Representation:</h4>
              <p>Angular Component: {result.tpExpression.angular}</p>
              <p>Prime Component: {result.tpExpression.prime}</p>
            </div>
            
            <div className="mt-4 border p-4 rounded bg-green-50">
              <h4 className="font-medium mb-2">Decoded Polynomial:</h4>
              <p>{result.decoded}</p>
            </div>
          </div>
        );
      
      case 'derivative':
        return (
          <div className="bg-white p-6 rounded-lg shadow-md">
            <h3 className="text-xl font-semibold mb-4">Derivative Result:</h3>
            <p className="mb-4">f(x) = {result.operand}</p>
            <p className="mb-4">f'(x) = {result.result}</p>
            
            <ResponsiveContainer width="100%" height={300}>
              <LineChart data={result.data} margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
                <CartesianGrid strokeDasharray="3 3" />
                <XAxis dataKey="x" />
                <YAxis />
                <Tooltip />
                <Legend />
                <Line type="monotone" dataKey="original" stroke="#8884d8" name="f(x)" dot={false} />
                <Line type="monotone" dataKey="derived" stroke="#82ca9d" name="f'(x)" dot={false} />
              </LineChart>
            </ResponsiveContainer>
            
            <div className="mt-6 grid grid-cols-1 md:grid-cols-2 gap-4">
              <div className="border p-4 rounded bg-blue-50">
                <h4 className="font-medium mb-2">TrigoPrime Representation of f(x):</h4>
                <p>Angular Component: {result.tpExpression.angular}</p>
                <p>Prime Component: {result.tpExpression.prime}</p>
              </div>
              
              <div className="border p-4 rounded bg-blue-50">
                <h4 className="font-medium mb-2">TrigoPrime Representation of f'(x):</h4>
                <p>Angular Component: {result.tpDerivative.angular}</p>
                <p>Prime Component: {result.tpDerivative.prime}</p>
              </div>
            </div>
            
            <div className="mt-4 border p-4 rounded bg-green-50">
              <h4 className="font-medium mb-2">Decoded Derivative:</h4>
              <p>{result.decoded}</p>
            </div>
          </div>
        );
      
      case 'differential-equation':
        return (
          <div className="bg-white p-6 rounded-lg shadow-md">
            <h3 className="text-xl font-semibold mb-4">Differential Equation Solution:</h3>
            <p className="mb-4">Equation: dy/dx = {result.operand}</p>
            <p className="mb-4">Solution: {result.result}</p>
            
            <ResponsiveContainer width="100%" height={300}>
              <LineChart data={result.numerical} margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
                <CartesianGrid strokeDasharray="3 3" />
                <XAxis dataKey="x" />
                <YAxis />
                <Tooltip />
                <Line type="monotone" dataKey="y" stroke="#8884d8" name="y(x)" dot={false} />
              </LineChart>
            </ResponsiveContainer>
            
            <div className="mt-6 grid grid-cols-1 md:grid-cols-2 gap-4">
              <div className="border p-4 rounded bg-blue-50">
                <h4 className="font-medium mb-2">TrigoPrime Representation of Equation:</h4>
                <p>Angular Component: {result.tpEquation.angular}</p>
                <p>Prime Component: {result.tpEquation.prime}</p>
              </div>
              
              <div className="border p-4 rounded bg-blue-50">
                <h4 className="font-medium mb-2">TrigoPrime Representation of Solution:</h4>
                <p>Angular Component: {result.tpSolution.angular}</p>
                <p>Prime Component: {result.tpSolution.prime}</p>
              </div>
            </div>
            
            <div className="mt-4 border p-4 rounded bg-green-50">
              <h4 className="font-medium mb-2">Decoded Solution:</h4>
              <p>{result.decoded}</p>
            </div>
          </div>
        );
      
      case 'fourier-series':
        return (
          <div className="bg-white p-6 rounded-lg shadow-md">
            <h3 className="text-xl font-semibold mb-4">Fourier Series Approximation:</h3>
            <p className="mb-4">Original Function: f(x) = {result.originalFunction}</p>
            <p className="mb-4">Fourier Series: {result.fourierSeries}</p>
            
            <FourierSeriesVisualizer 
              originalFn={result.originalFunction} 
              terms={result.terms} 
              domain={[-Math.PI, Math.PI]} 
              title="Fourier Series Approximation" 
            />
            
            <div className="mt-6 grid grid-cols-1 md:grid-cols-2 gap-4">
              <div className="border p-4 rounded bg-blue-50">
                <h4 className="font-medium mb-2">TrigoPrime Representation of Original Function:</h4>
                <p>Angular Component: {result.tpOriginalFunction.angular}</p>
                <p>Prime Component: {result.tpOriginalFunction.prime}</p>
              </div>
              
              <div className="border p-4 rounded bg-blue-50">
                <h4 className="font-medium mb-2">TrigoPrime Representation of Fourier Series:</h4>
                <p>Angular Component: {result.tpFourierSeries.angular}</p>
                <p>Prime Component: {result.tpFourierSeries.prime}</p>
              </div>
            </div>
            
            <div className="mt-4 border p-4 rounded bg-green-50">
              <h4 className="font-medium mb-2">Decoded Fourier Series:</h4>
              <p>{result.decoded}</p>
            </div>
          </div>
        );
      
      case 'wave-function-analysis':
        return (
          <div className="bg-white p-6 rounded-lg shadow-md">
            <h3 className="text-xl font-semibold mb-4">Quantum Wave Function Visualization:</h3>
            <p className="mb-4">ψ(x) = {result.waveFunction}</p>
            
            <WaveFunctionVisualizer 
              psiRe={result.psiRe} 
              psiIm={result.psiIm} 
              domain={[-10, 10]} 
              title="Wave Function" 
            />
            
            <div className="mt-4 grid grid-cols-1 md:grid-cols-2 gap-4">
              <div className="border p-4 rounded bg-gray-50">
                <h4 className="font-medium mb-2">Wave Function Properties:</h4>
                <p>Maximum Probability: {result.maxProbability.toFixed(4)} at x = {result.maxProbabilityPosition.toFixed(4)}</p>
                <p>Total Probability: {result.totalProbability.toFixed(4)}</p>
                <p>Normalization Factor: {result.normalizationFactor.toFixed(4)}</p>
              </div>
              
              <div className="border p-4 rounded bg-blue-50">
                <h4 className="font-medium mb-2">TrigoPrime Representation:</h4>
                <p>Angular Component: {result.tpWaveFunction.angular}</p>
                <p>Prime Component: {result.tpWaveFunction.prime}</p>
              </div>
            </div>
            
            <div className="mt-4 border p-4 rounded bg-green-50">
              <h4 className="font-medium mb-2">Decoded Wave Function:</h4>
              <p>{result.decoded}</p>
            </div>
          </div>
        );
      
      default:
        return <div className="text-red-500 p-4">Unknown operation type</div>;
    }
  };
  
  return (
    <div className="min-h-screen bg-gray-100 p-4">
      <div className="max-w-6xl mx-auto">
        <h1 className="text-3xl font-bold text-center mb-8">TrigoPrime Symbolic Solver</h1>
        
        {/* Navigation Tabs */}
        <div className="flex flex-wrap mb-6">
          <button 
            className={`px-4 py-2 mr-2 mb-2 rounded-md ${activeTab === 'arithmetic' ? 'bg-blue-600 text-white' : 'bg-white'}`}
            onClick={() => setActiveTab('arithmetic')}
          >
            Arithmetic
          </button>
          <button 
            className={`px-4 py-2 mr-2 mb-2 rounded-md ${activeTab === 'matrices' ? 'bg-blue-600 text-white' : 'bg-white'}`}
            onClick={() => setActiveTab('matrices')}
          >
            Matrices
          </button>
          <button 
            className={`px-4 py-2 mr-2 mb-2 rounded-md ${activeTab === 'polynomials' ? 'bg-blue-600 text-white' : 'bg-white'}`}
            onClick={() => setActiveTab('polynomials')}
          >
            Polynomials
          </button>
          <button 
            className={`px-4 py-2 mr-2 mb-2 rounded-md ${activeTab === 'calculus' ? 'bg-blue-600 text-white' : 'bg-white'}`}
            onClick={() => setActiveTab('calculus')}
          >
            Calculus
          </button>
          <button 
            className={`px-4 py-2 mr-2 mb-2 rounded-md ${activeTab === 'diffeq' ? 'bg-blue-600 text-white' : 'bg-white'}`}
            onClick={() => setActiveTab('diffeq')}
          >
            Differential Equations
          </button>
          <button 
            className={`px-4 py-2 mr-2 mb-2 rounded-md ${activeTab === 'fourier' ? 'bg-blue-600 text-white' : 'bg-white'}`}
            onClick={() => setActiveTab('fourier')}
          >
            Fourier Analysis
          </button>
          <button 
            className={`px-4 py-2 mr-2 mb-2 rounded-md ${activeTab === 'quantum' ? 'bg-blue-600 text-white' : 'bg-white'}`}
            onClick={() => setActiveTab('quantum')}
          >
            Quantum Mechanics
          </button>
        </div>
        
        {/* Input Forms */}
        <div className="bg-white p-6 rounded-lg shadow-md mb-6">
          {/* Arithmetic Tab */}
          {activeTab === 'arithmetic' && (
            <div>
              <h2 className="text-2xl font-semibold mb-4">Basic Arithmetic Operations</h2>
              <div className="mb-4">
                <label className="block mb-2">Value A:</label>
                <input 
                  type="text" 
                  value={inputA} 
                  onChange={(e) => setInputA(e.target.value)} 
                  className="w-full p-2 border rounded-md"
                />
              </div>
              <div className="mb-4">
                <label className="block mb-2">Value B:</label>
                <input 
                  type="text" 
                  value={inputB} 
                  onChange={(e) => setInputB(e.target.value)} 
                  className="w-full p-2 border rounded-md"
                />
              </div>
              <div className="flex space-x-4">
                <button 
                  className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700" 
                  onClick={() => handleArithmeticSubmit('add')}
                >
                  Add (A + B)
                </button>
                <button 
                  className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700" 
                  onClick={() => handleArithmeticSubmit('multiply')}
                >
                  Multiply (A × B)
                </button>
              </div>
            </div>
          )}
          
          {/* Matrices Tab */}
          {activeTab === 'matrices' && (
            <div>
              <h2 className="text-2xl font-semibold mb-4">Matrix Operations</h2>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mb-4">
                <div>
                  <h3 className="text-lg font-medium mb-2">Matrix A:</h3>
                  <div className="flex flex-col items-center">
                    {matrixA.map((row, i) => (
                      <div key={`row-a-${i}`} className="flex mb-2">
                        {row.map((cell, j) => (
                          <input
                            key={`cell-a-${i}-${j}`}
                            type="number"
                            value={cell}
                            onChange={(e) => handleMatrixCellChange(matrixA, setMatrixA, i, j, e.target.value)}
                            className="w-16 h-16 text-center border m-1"
                          />
                        ))}
                      </div>
                    ))}
                  </div>
                </div>
                <div>
                  <h3 className="text-lg font-medium mb-2">Matrix B:</h3>
                  <div className="flex flex-col items-center">
                    {matrixB.map((row, i) => (
                      <div key={`row-b-${i}`} className="flex mb-2">
                        {row.map((cell, j) => (
                          <input
                            key={`cell-b-${i}-${j}`}
                            type="number"
                            value={cell}
                            onChange={(e) => handleMatrixCellChange(matrixB, setMatrixB, i, j, e.target.value)}
                            className="w-16 h-16 text-center border m-1"
                          />
                        ))}
                      </div>
                    ))}
                  </div>
                </div>
              </div>
              <div className="flex space-x-4">
                <button 
                  className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700" 
                  onClick={() => handleMatrixSubmit('add')}
                >
                  Add (A + B)
                </button>
                <button 
                  className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700" 
                  onClick={() => handleMatrixSubmit('multiply')}
                >
                  Multiply (A × B)
                </button>
              </div>
            </div>
          )}
          
          {/* Polynomials Tab */}
          {activeTab === 'polynomials' && (
            <div>
              <h2 className="text-2xl font-semibold mb-4">Polynomial Operations</h2>
              <div className="mb-4">
                <label className="block mb-2">Polynomial Expression (in terms of x):</label>
                <input 
                  type="text" 
                  value={polynomialExpr} 
                  onChange={(e) => setPolynomialExpr(e.target.value)} 
                  className="w-full p-2 border rounded-md"
                  placeholder="e.g., 2*x^2 + 3*x + 1"
                />
              </div>
              <button 
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700" 
                onClick={handlePolynomialSubmit}
              >
                Visualize Polynomial
              </button>
            </div>
          )}
          
          {/* Calculus Tab */}
          {activeTab === 'calculus' && (
            <div>
              <h2 className="text-2xl font-semibold mb-4">Calculus Operations</h2>
              <div className="mb-4">
                <label className="block mb-2">Function Expression (in terms of x):</label>
                <input 
                  type="text" 
                  value={calcExpr} 
                  onChange={(e) => setCalcExpr(e.target.value)} 
                  className="w-full p-2 border rounded-md"
                  placeholder="e.g., x^2"
                />
              </div>
              <button 
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700" 
                onClick={handleCalculusSubmit}
              >
                Compute Derivative
              </button>
            </div>
          )}
          
          {/* Differential Equations Tab */}
          {activeTab === 'diffeq' && (
            <div>
              <h2 className="text-2xl font-semibold mb-4">Differential Equations</h2>
              <div className="mb-4">
                <label className="block mb-2">Differential Equation (dy/dx = ...):</label>
                <input 
                  type="text" 
                  value={diffEqExpr} 
                  onChange={(e) => setDiffEqExpr(e.target.value)} 
                  className="w-full p-2 border rounded-md"
                  placeholder="e.g., y"
                />
              </div>
              <button 
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700" 
                onClick={handleDiffEqSubmit}
              >
                Solve Differential Equation
              </button>
            </div>
          )}
          
          {/* Fourier Analysis Tab */}
          {activeTab === 'fourier' && (
            <div>
              <h2 className="text-2xl font-semibold mb-4">Fourier Series Analysis</h2>
              <div className="mb-4">
                <label className="block mb-2">Function Expression (in terms of x):</label>
                <input 
                  type="text" 
                  value={fourierExpr} 
                  onChange={(e) => setFourierExpr(e.target.value)} 
                  className="w-full p-2 border rounded-md"
                  placeholder="e.g., x"
                />
              </div>
              <div className="mb-4">
                <label className="block mb-2">Number of Terms:</label>
                <input 
                  type="number" 
                  value={fourierTerms} 
                  onChange={(e) => setFourierTerms(parseInt(e.target.value) || 1)} 
                  className="w-full p-2 border rounded-md"
                  min="1"
                  max="50"
                />
              </div>
              <button 
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700" 
                onClick={handleFourierSubmit}
              >
                Compute Fourier Series
              </button>
            </div>
          )}
          
          {/* Quantum Mechanics Tab */}
          {activeTab === 'quantum' && (
            <div>
              <h2 className="text-2xl font-semibold mb-4">Quantum Wave Functions</h2>
              <div className="mb-4">
                <label className="block mb-2">Real Part Expression (in terms of x):</label>
                <input 
                  type="text" 
                  value={waveReExpr} 
                  onChange={(e) => setWaveReExpr(e.target.value)} 
                  className="w-full p-2 border rounded-md"
                  placeholder="e.g., cos(x)"
                />
              </div>
              <div className="mb-4">
                <label className="block mb-2">Imaginary Part Expression (in terms of x):</label>
                <input 
                  type="text" 
                  value={waveImExpr} 
                  onChange={(e) => setWaveImExpr(e.target.value)} 
                  className="w-full p-2 border rounded-md"
                  placeholder="e.g., sin(x)"
                />
              </div>
              <button 
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700" 
                onClick={handleWaveFunctionSubmit}
              >
                Analyze Wave Function
              </button>
            </div>
          )}
        </div>
        
        {/* Results Section */}
        {result && (
          <div className="mb-6">
            {renderResult()}
          </div>
        )}
        
        {/* Documentation Section */}
        <div className="bg-white p-6 rounded-lg shadow-md">
          <h2 className="text-2xl font-semibold mb-4">TrigoPrime Framework</h2>
          <p className="mb-4">
            The TrigoPrime framework is a dual-representation system that combines complex phase representation (Angular Set Theory) 
            and prime factorization encoding to create a unified computational approach across mathematical domains.
          </p>
          <p className="mb-4">
            <strong>Core Definition:</strong> TP(x) = {"{e^(iπx), E(x)}"} where e^(iπx) represents the angular encoding and E(x) represents the prime encoding.
          </p>
          <p className="mb-4">
            <strong>Key Properties:</strong>
          </p>
          <ul className="list-disc ml-6 mb-4">
            <li>Bijective Mapping: TP(a) = TP(b) ⟺ a = b</li>
            <li>Addition: TP(a + b) = TP(a) ⊕ TP(b) = {"{e^(iπ(a+b)), E(a+b)}"}</li>
            <li>Multiplication: TP(a·b) = TP(a) ⊗ TP(b) = {"{e^(iπa)·e^(iπb), E(a)·E(b)}"}</li>
          </ul>
          <p>
            This solver allows you to explore these operations across various mathematical domains, visualize the results,
            and observe the transformations between original mathematical objects and their TrigoPrime representations.
          </p>
        </div>
      </div>
    </div>
  );
}

export default TrigoprimeApp;
