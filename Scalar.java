/**
 * Numeric solver for scalar ODE's
 *          Forward Euler's
 *          Huen's
 *          Runge-Kutta4
 *          
 *          
 * @author Trenton Sorenen
 *
 */
public class Scalar {
    /**
     * a set of numeric algorithms to be used on scaler differential equations
     * 
     * @param args
     *            not needed for this purpose
     */
    public static void main(String args[]) {
        System.out.println("Forward Euler");
        System.out.println("y(0.5) is approximately: " + forwardEuler(0, 0.5, 0.125,
            0.5));
        System.out.println("Heun's Method");
        System.out.println("y(0.5) is approximately: " + heunsMethod(0, 0.5, 0.25,
            0.5));
        System.out.println("4th order Runge-Kutta");
        System.out.println("y(0.5) is approximately: " + rungeKutta4(0, 0.5, 0.5));

    }


    /**
     * A method to carry out the forward Euler method
     * 
     * @param t0
     *            the initial time
     * @param y0
     *            the initial y position
     * @param h
     *            the time step
     * @param t
     *            the target time
     * @return an approximate solution to the differential equation
     */
    public static double forwardEuler(
        double t0,
        double y0,
        double h,
        double t) {
        double result = 0;
        // loops through the method the number of times that it takes the
        // initial
        // time plus the time step to reach the target time, t0 is updated each
        // iteration
        for (int i = 0; t0 < t; i++) {
            result = y0 + h * f(y0, t0);
            t0 = t0 + h;
            y0 = result;
        }

        return result;
    }


    /**
     * A method to carry out Heune's Method
     * 
     * @param t0
     *            the initial time
     * @param y0
     *            the initial y position
     * @param h
     *            the time step
     * @param t
     *            the target time
     * @return an approximate solution to the differential equation
     */
    public static double heunsMethod(double t0, double y0, double h, double t) {
        double result = 0;
        // loops through the method the number of times that it takes the
        // initial
        // time plus the time step to reach the target time, t0 is updated each
        // iteration
        for (int i = 0; t0 < t; i++) {
            double t1 = t0 + h;
            result = y0 + (h / 2) * ((f(y0, t0) + f(y0 + h * f(y0, t0), t1)));
            t0 = t1;
            y0 = result;
        }

        return result;

    }


    /**
     * A method to carry out Runge-Kutta4 method. A target time is not used in this implementation because h = t
     * @param t0
     *            the initial time
     * @param y0
     *            the initial y position
     * @param h
     *            the time step
     * @return an approximate solution to the differential equation
     */
    public static double rungeKutta4(double t0, double y0, double h) {
        double k1 = f(y0, t0);
        double k2 = f(y0 + ((h * k1) / 2), t0 + (h / 2));
        double k3 = f(y0 + ((h * k2) / 2), t0 + (h / 2));
        double k4 = f(y0 + (h * k3), t0 + h);

        return y0 + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

    }

/**
 * a method to set up and use the function that is in the differential equation.
 * @param y the y position
 * @param t the time 
 * @return the function that is being evaluated
 */
    public static double f(double y, double t) {
        return 4 * t * y * y;
    }
}
