/**
 * Numeric solver for scalar ODE's
 *          Forward Euler's
 *          Backwards Euler's
 *          Huen's
 *          Runge-Kutta4
 *          
 *          
 * @author Trenton Sorenen
 *
 */
public class Linear {
    /**
     * a set of numeric methods that can be used for linear differential
     * equations
     * 
     * @param args
     *            not used here
     */
    public static void main(String args[]) {

        double[] y = new double[2];
        y[0] = 1;
        y[1] = 1;
        double t0 = 0;
        System.out.println("Forward Euler");
        System.out.println("y(0.1) is approximately " + forwardEuler(y, t0, 0.1,
            0.1)[0] + ", " + forwardEuler(y, t0, 0.1, 0.1)[1]);
        System.out.println("Heun's Method");
        System.out.println("y(0.1) is approximately " + heunsMethod(y, t0, 0.1,
            0.1)[0] + ", " + heunsMethod(y, t0, 0.25, 0.1)[1]);
        System.out.println("Backwordward Euler");
        System.out.println("y(0.1) is approximately " + backwardEuler(y, t0, 0.1,
            0.1)[0] + ", " + backwardEuler(y, t0, 0.1, 0.1)[1]);
        System.out.println("4th order Runge-Kutta");
        System.out.println("y(0.1) is approximately " + rungeKutta4(y, t0, 0.1)[0]
            + ", " + rungeKutta4(y, t0, 0.1)[1]);

    }


    /**
     * A method to carry out Forward Euler on a linear differential equation
     * 
     * @param y0
     *            the initial y positions
     * @param t0
     *            the initial time
     * @param h
     *            the time step
     * @param t
     *            the target time
     * @return an approximation to the differential equation.
     */
    public static double[] forwardEuler(
        double[] y0,
        double t0,
        double h,
        double t) {
        int l = y0.length;
        // loops through the method the number of time corresponding to the
        // number of elements in y0
        double[] result = new double[2];
        for (int i = 0; i < l - 1; i++) {
            result[i] = y0[i] + h * f(y0, t0)[i];
            result[i + 1] = y0[i + 1] + h * f(y0, t0)[i + 1];
            t0 = t0 + h;
            y0 = result;
        }

        return result;

    }


    /**
     * a method to carry out the Backward Euler method on a linear differential equation
     * @param y0
     *            the initial y positions
     * @param t0
     *            the initial time
     * @param h
     *            the time step
     * @param t
     *            the target time
     * @return an approximation to the differential equation.
     */
    public static double[] backwardEuler(
        double[] y0,
        double t0,
        double h,
        double t) {
        //the A matrix of f(y,t)
        int l = y0.length;
        double[][] A = new double[2][2];
        A[0][0] = -2;
        A[0][1] = t;
        A[1][0] = 3;
        A[1][1] = 0;
        // an Identity matrix
        double[][] I = new double[2][2];
        I[0][0] = 1;
        I[1][1] = 1;
        I[0][1] = 0;
        I[1][0] = 0;
        double[][] C = new double[2][2];
        // carries out (I-h*A)
        for (int i = 0; i < l; i++) {
            for (int j = 0; j < l; j++)
                C[i][j] = (I[i][j] - h * A[i][j]);
        }
        // solve the 2x2 matrix system C*y1=y0
        double[] y = new double[l];
        y[1] = y0[1] / (((C[1][0] * y0[0] - C[0][1]) / C[0][0]) + C[1][1]);
        y[0] = (y0[0] - C[0][1] * y[1]) / C[0][0];

        return y;

    }

    
    /**
     * a method to carry out Heune's method on a linear differential equation.
     * This was easiest done as an RK-2 implementation to keep track of the
     * functions with in functions. yk+1= y_k +h/2(f(t_k,y_k)+f(t_(k+1),y_k+h*f(t_k,y_k)))
     * @param y0
     *            the initial y positions
     * @param t0
     *            the initial time
     * @param h
     *            the time step
     * @param t
     *            the target time
     * @return an approximation to the differential equation.
     */

    public static double[] heunsMethod(
        double[] y0,
        double t0,
        double h,
        double t) {

        int l = y0.length;
        double[] k1 = new double[l];
        double[] k2 = new double[l];
        k1 = f(y0, t0);
        for (int i = 0; i < l; i++) {
            k2[i] = y0[i] + h * k1[i] / 2;
        }
        k2 = f(k2, t + h / 2);

        for (int i = 0; i < l; i++) {

            k1[i] = y0[i] + (h / 2 * (k1[i])) + (h / 2 * (k2[i]));
        }
        return k1;

    }

    /**
     * a method to carry out Runge-Kutta 4 on a linear differential equation.
     * 
     * @param y0
     *            the initial y positions
     * @param t0
     *            the initial time
     * @param h
     *            the time step
     * @param t
     *            the target time
     * @return an approximation to the differential equation.
     */
    public static double[] rungeKutta4(double y0[], double t, double h) {
        int l = y0.length;
        double k1[] = new double[l];
        double k2[] = new double[l];
        double k3[] = new double[l];
        double k4[] = new double[l];
        k1 = f(y0, t);
        for (int i = 0; i < l; ++i) {
            k2[i] = y0[i] + h * k1[i] / 2;
        }
        k2 = f(k2, t + h / 2);

        for (int i = 0; i < l; ++i) {
            k3[i] = y0[i] + h * k2[i] / 2;
        }
        k3 = f(k3, t + h / 2);
        for (int i = 0; i < l; ++i) {
            k4[i] = y0[i] + h * k3[i];
        }
        k4 = f(k4, t + h);
        for (int i = 0; i < l; ++i) {
            k1[i] = y0[i] + h * (k1[i] + 2 * (k2[i] + k3[i]) + k4[i]) / 6;
        }
        return k1;
    }


    /**
     * sets up the linear function
     * 
     * @param y
     * @param t
     * @return
     */
    public static double[] f(double[] y, double t) {
        double[][] A = new double[2][2];
        A[0][0] = -2;
        A[0][1] = t;
        A[1][0] = 3;
        A[1][1] = 0;
        // double[] Y = new double[2];
// Y[1] = y;
// Y[2] = y;
        return matrixMultiply(A, y);
    }


    /**
     * multiplies a 2x2 * 2x1
     * 
     * @param A
     *            a 2x2 matrix
     * @param y
     *            a 2x1 vector
     * @return
     */
    public static double[] matrixMultiply(double[][] A, double[] y) {
        double[] result = new double[2];
        result[0] = A[0][0] * y[0] + A[0][1] * y[1];
        result[1] = A[1][0] * y[0] + A[1][1] * y[1];

        return result;

    }
}
