public class Matrix {

    private final int height, width;
    private final double[][] matrix;

    public Matrix(double[][] matrix) {
        height = matrix.length;
        width = matrix[0].length;
        this.matrix = matrix;
    }
    
    public Matrix(int height, int width) {
        this(new double[height][width]);
    }
    
    public double[][] getArray() {
        return matrix;
    }
    
     public double[][] getArrayCopy () {
      double[][] C = new double[height][width];
      for (int i = 0; i < height; i++) {
         for (int j = 0; j < width; j++) {
            C[i][j] = matrix[i][j];
         }
      }
      return C;
    }


    public double get(int i, int j) {
        return matrix[i][j];
    }

    public int getHeight() {
        return height;
    }

    public int getWidth() {
        return width;
    }
    
    public double[] getRow(int i1) {
        double[] toreturn = new double[width];
        for (int i = 0; i < width; i++) {
            toreturn[i] = matrix[i1][i];
        }
        return toreturn;
    }
    
    public double[] getColumn(int i1) {
        double[] toreturn = new double[height];
        for (int i = 0; i < height; i++) {
            toreturn[i] = matrix[i][i1];
        }
        return toreturn;
    }
    
    public Matrix getMatrix (int i0, int i1, int j0, int j1) {
      Matrix X = new Matrix(i1-i0+1,j1-j0+1);
      double[][] B = X.getArray();
      try {
         for (int i = i0; i <= i1; i++) {
            for (int j = j0; j <= j1; j++) {
               B[i-i0][j-j0] = matrix[i][j];
            }
         }
      } catch(ArrayIndexOutOfBoundsException e) {
         throw new ArrayIndexOutOfBoundsException("Submatrix indices");
      }
      return X;
   }
    
    public Matrix getMatrix(int[] r, int cInit, int cFinal) { 
        int cNum = cFinal - cInit + 1;
        double[][] temp = new double[r.length][cNum];
        Matrix x = new Matrix(temp);
        double[][] b = x.getArray();
        
        try {
            for (int i = 0; i <= r.length; i++) {
                for (int j = cInit; j <= cFinal; j++) {
                    b[i][j-cInit] = matrix[i][j];
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return x;
    }
    
    public void setValue(int row, int col, double val) {
        matrix[row][col] = val;
    }
    
    public Matrix transpose() {
        double transposed[][] = new double[width][height];
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                transposed[j][i] = matrix[i][j];
            }
        }
        Matrix toreturn = new Matrix(transposed);
        return toreturn;
    }

    public Matrix inverse2x2() {
        double adj[][] = new double[width][height];
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                adj[i][j] = matrix[i][j];
            }
        }
        double A = adj[0][0];
        double D = adj[1][1];
        adj[0][0] = D;
        adj[1][1] = A;
        adj[0][1] = -(adj[0][1]);
        adj[1][0] = -(adj[1][0]);

        double detOfA = LinearAlgebra.determinant(new Matrix(matrix));
        Matrix adjMatrix = new Matrix(adj);
        Matrix inverse = LinearAlgebra.matrixScalarMultiple(adjMatrix, (1 / detOfA));
        return inverse;
    }

    public String toString() {
        String s = new String("");
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                s += matrix[r][c] + "   ";
            }
            s += "\n";
        }
        return s;
    }

    public double trace() {
        double temp[][] = new double[width][height];
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                temp[i][j] = matrix[i][j];
            }
        }
        double trace = 0;
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                if (i == j) {
                     trace = trace + temp[i][j];
                }
            }
        }
        return trace;
    }
}