package sph;

public class Particle {
	double mass = 0.1;
	double rho, e, ePred, c, eDot, P ;
	double[] r;
	double [] v;
	double [] vPred;
	double [] a;
	double vAbs;
	
	public PrioQ pq;
	
	public Particle(double mass, double x, double y) {
		pq = new PrioQ();
		r = new double[2];
		this.mass = mass;
		r[0] = x;
		r[1] = y;
		v = new double[2];
		v[0] = 0; v[1] = 0;
		vPred = new double[2];
		a = new double[2];
		e = 2;
	}
}
