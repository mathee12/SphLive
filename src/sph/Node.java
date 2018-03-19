package sph;

public class Node {
	Node left, right;
	public int start, end;
	public double[] rmin, rmax;
	private boolean leftTrue = true;
	public double balancePoint[];
	public int nodeParticles;

	public Node(int start, int end, double xmin, double xmax, double ymin, double ymax) {
		rmin = new double[2];
		rmax = new double[2];
		rmin[0] = xmin;
		rmin[1] = ymin;
		rmax[0] = xmax;
		rmax[1] = ymax;
		this.start = start;
		this.end = end;
		balancePoint = new double[2];
		nodeParticles = end - start + 1;
	}
	
	public void setLeftTrue(boolean bool) {
		this.leftTrue = bool;
	}
	
	public boolean getLeftTrue() {
		return leftTrue;
	}
}
