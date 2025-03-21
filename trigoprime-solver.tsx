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
                    backgroundColor: `rgba(66, 153, 225, ${Math.abs(cell/10)})`,
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

// Component for visualizing Fourier Series
const FourierSeriesVisualizer = ({ originalFn, terms = 5, domain = [-Math.PI, Math.PI], points = 200, title }) => {
  const [seriesData, setSeriesData] = useState([]);
  
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
  // Basic conversion functions
  static toTrigoRepresentation(x) {
    if (typeof x === 'number') {
      // For integers
      if (Number.isInteger(x)) {
        const angular = Math.exp(Math.PI * x * Math.I);
        const prime = this.getPrimeEncoding(x);
        return { angular, prime, original: x };
      }
      // For rational numbers
      else if (isFinite(x)) {
        const [num, den] = this.toFraction(x);
        const angular = Math.exp(Math.PI * (num/den) * Math.I);
        const prime = this.getPrimeEncoding(num) / this.getPrimeEncoding(den);
        return { angular, prime, original: x };
      }
    }
    // For complex numbers
    else if (x.re !== undefined && x.im !== undefined) {
      const magnitude = Math.sqrt(x.re**2 + x.im**2);
      const angle = Math.atan2(x.im, x.re);
      const angular = Math.exp(Math.PI * magnitude * Math.exp(angle * Math.I) * Math.I);
      // Simplified prime encoding for complex numbers
      const prime = `${this.getPrimeEncoding(x.re)} * E(i)^${this.getPrimeEncoding(x.im)}`;
      return { angular, prime, original: x };
    }
    
    return { angular: null, prime: null, original: x };
  }
  
  // Helper function to get prime encoding (simplified for demonstration)
  static getPrimeEncoding(n) {
    if (n === 0) return 1;
    
    // Get the nth prime for positive integers
    if (n > 0 && Number.isInteger(n)) {
      const primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
      return n <= primes.length ? primes[n-1] : `Prime(${n})`;
    }
    
    // For negative integers
    if (n < 0 && Number.isInteger(n)) {
      return 2 * this.getPrimeEncoding(-n);
    }
    
    // For other numbers, return a representation
    return `E(${n})`;
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
    } while (Math.abs(decimal - h1 / k1) > decimal * tolerance);
    
    return [h1, k1];
  }
  
  // Addition in TrigoPrime
  static add(a, b) {
    const aTP = this.toTrigoRepresentation(a);
    const bTP = this.toTrigoRepresentation(b);
    
    // Direct addition for the original values
    const result = aTP.original + bTP.original;
    
    // TP representation of the result
    const resultTP = this.toTrigoRepresentation(result);
    
    return {
      result,
      tpResult: resultTP,
      operation: 'addition',
      operands: [aTP, bTP]
    };
  }
  
  // Multiplication in TrigoPrime
  static multiply(a, b) {
    const aTP = this.toTrigoRepresentation(a);
    const bTP = this.toTrigoRepresentation(b);
    
    // Direct multiplication for the original values
    const result = aTP.original * bTP.original;
    
    // TP representation of the result
    const resultTP = this.toTrigoRepresentation(result);
    
    return {
      result,
      tpResult: resultTP,
      operation: 'multiplication',
      operands: [aTP, bTP]
    };
  }
  
  // Matrix operations
  static matrixAdd(A, B) {
    if (!A || !B || A.length !== B.length || A[0].length !== B[0].length) {
      throw new Error("Matrices must have the same dimensions for addition");
    }
    
    const result = [];
    for (let i = 0; i < A.length; i++) {
      result[i] = [];
      for (let j = 0; j < A[0].length; j++) {
        result[i][j] = A[i][j] + B[i][j];
      }
    }
    
    // Convert to TrigoRepresentation
    const tpResult = result.map(row => row.map(val => this.toTrigoRepresentation(val)));
    
    return {
      result,
      tpResult,
      operation: 'matrix-addition',
      operands: [A, B]
    };
  }
  
  static matrixMultiply(A, B) {
    if (!A || !B || A[0].length !== B.length) {
      throw new Error("Matrix dimensions are not compatible for multiplication");
    }
    
    const result = [];
    for (let i = 0; i < A.length; i++) {
      result[i] = [];
      for (let j = 0; j < B[0].length; j++) {
        result[i][j] = 0;
        for (let k = 0; k < A[0].length; k++) {
          result[i][j] += A[i][k] * B[k][j];
        }
      }
    }
    
    // Convert to TrigoRepresentation
    const tpResult = result.map(row => row.map(val => this.toTrigoRepresentation(val)));
    
    return {
      result,
      tpResult,
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
      
      return {
        result: derivative,
        data,
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
    return {
      result: "Analytical solution placeholder",
      numerical: this.generateNumericalSolution(equation, initialConditions),
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
      
      // For demonstration, we'll use a simple ODE: dy/dx = y (solution: y = e^x)
      const slope = y; // Replace with actual equation evaluation
      y = y + h * slope;
      x = x + h;
    }
    
    return data;
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
        expression: polynomialExpr
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
      setResult({
        operation: 'fourier',
        expression: fourierExpr,
        terms: fourierTerms
      });
    } catch (error) {
      console.error("Error:", error);
      setResult({ error: error.message });
    }
  };
  
  const handleWaveFunctionSubmit = () => {
    try {
      setResult({
        operation: 'wave-function',
        reExpr: waveReExpr,
        imExpr: waveImExpr
      });
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
            <div className="grid grid-cols-2 gap-4">
              <div>
                <h4 className="font-medium mb-2">Original Values:</h4>
                <p>a = {result.operands[0].original}</p>
                <p>b = {result.operands[1].original}</p>
              </div>
              <div>
                <h4 className="font-medium mb-2">TrigoPrime Representation:</h4>
                <p>TP(result) = {`{${result.tpResult.angular}, ${result.tpResult.prime}}`}</p>
              </div>
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
          </div>
        );
      
      case 'polynomial':
        return (
          <div className="bg-white p-6 rounded-lg shadow-md">
            <h3 className="text-xl font-semibold mb-4">Polynomial Visualization:</h3>
            <PolynomialPlot polynomial={result.expression} domain={[-5, 5]} title={`P(x) = ${result.expression}`} />
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
          </div>
        );
      
      case 'fourier':
        return (
          <div className="bg-white p-6 rounded-lg shadow-md">
            <h3 className="text-xl font-semibold mb-4">Fourier Series Approximation:</h3>
            <p className="mb-4">Original Function: f(x) = {result.expression}</p>
            <p className="mb-4">Number of Terms: {result.terms}</p>
            <FourierSeriesVisualizer 
              originalFn={result.expression} 
              terms={result.terms} 
              domain={[-Math.PI, Math.PI]} 
              title="Fourier Series Approximation" 
            />
          </div>
        );
      
      case 'wave-function':
        return (
          <div className="bg-white p-6 rounded-lg shadow-md">
            <h3 className="text-xl font-semibold mb-4">Quantum Wave Function Visualization:</h3>
            <p className="mb-4">Re(ψ) = {result.reExpr}</p>
            <p className="mb-4">Im(ψ) = {result.imExpr}</p>
            <WaveFunctionVisualizer 
              psiRe={result.reExpr} 
              psiIm={result.imExpr} 
              domain={[-10, 10]} 
              title="Wave Function" 
            />
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
                Visualize Wave Function
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
            Core Definition: TP(x) = {"{e^(iπx), E(x)}"} where e^(iπx) represents the angular encoding and E(x) represents the prime encoding.
          </p>
          <p>
            This solver allows you to explore various mathematical operations within the TrigoPrime framework and visualize the results.
          </p>
        </div>
      </div>
    </div>
  );
}

export default TrigoprimeApp;
