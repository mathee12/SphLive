package sph;

public class PrioQ {
	int index;
	double[] dist2, listdx, listdy;
	int[] list;
	public int NN = 32;
	
	public void init() {
		dist2 = new double[NN];
		list = new int[NN];
		listdx = new double[NN];
		listdy = new double[NN];
		for (int i = 0; i < NN; i++) {
			list[i] = SphLive.numberOfParticles;
			dist2[i] = Double.POSITIVE_INFINITY;
		}
		index = 0;
	}
	
	public double mindist2 (Node cell, Particle particle, double[] offset) {
		double dx1 = particle.r[0] - (cell.rmax[0] + offset[0]);
		double dx2 = cell.rmin[0] + offset[0] - particle.r[0];
		if (dx2 > dx1) {
			dx1 = dx2;
		}
		if (dx1 < 0) {
			dx1 = 0;
		}
		double dy1 = particle.r[1] - (cell.rmax[1] + offset[1]);
		double dy2 = cell.rmin[1] + offset[1] - particle.r[1];
		if (dy2 > dy1) {
			dy1 = dy2;
		}
		if (dy1 < 0) {
			dy1 = 0;
		}
		return dx1*dx1 + dy1*dy1;
	}
	
	public double mindist2 (Node cell, Particle particle) {
		double dx1 = particle.r[0] - cell.rmax[0];
		double dx2 = cell.rmin[0] - particle.r[0];
		if (dx2 > dx1) {
			dx1 = dx2;
		}
		if (dx1 < 0) {
			dx1 = 0;
		}
		double dy1 = particle.r[1] - cell.rmax[1];
		double dy2 = cell.rmin[1] - particle.r[1];
		if (dy2 > dy1) {
			dy1 = dy2;
		}
		if (dy1 < 0) {
			dy1 = 0;
		}
		return dx1*dx1 + dy1*dy1;
	}
	
	public void replace(double dx,double dy, int j) {
		double d2 = dx*dx + dy*dy;
		if (d2 < dist2[index]) {
			list[index] = j;
			dist2[index] = d2;
			listdx[index] = dx;			listdy[index] = dy;
		
			double max = 0;
			for (int i = 0; i < list.length; i++) {
				if (dist2[i] > max) {
					max = dist2[i];
					index = i;
				}
			}
		}
	}
	
	public double getRadius() {
		return Math.sqrt(dist2[index]);
	}
}
