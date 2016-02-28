public class Vector {

    private double[] vector;
    private final int length;
    private boolean isVertical;

    public Vector(double[] vector) {
        this.vector = vector;
        length = vector.length;
        isVertical = true;
    }
    
    public Vector(int size) {
        double[] v = new double[size];
        this.vector = v;
        length = size;
        isVertical = true;
    }
    
    public void setValue(int loc, double val) {
        vector[loc] = val;
    }
    
    public double norm() {
        double norm = 0;
        for (int i = 0; i < vector.length; i++) {
            norm += vector[i] * vector[i];
        }
        norm = Math.sqrt(norm);
        return norm;
    }

    public double get(int i) {
        return vector[i];
    }
    
    public boolean getIsVertical() {
        return isVertical;
    }
    
    public void setIsVertical(boolean b) {
        isVertical = b;
    }

    public int getLength() {
        return length;
    }

    public String toString() {
        String s = new String("");
        for (int i = 0; i < length; i++) {
            s += vector[i] + "  ";
        }
        return s;
    }
    
    public Vector transpose() {
        setIsVertical(!isVertical);
        return this;
    }
    
    public void setRow(double[] newRow) {
        vector = newRow;
    } 
}