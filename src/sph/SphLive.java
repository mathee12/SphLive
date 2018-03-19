package sph;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Toolkit;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;


@SuppressWarnings("serial")
public class SphLive extends JPanel{

	// ---------- PARAMETER SECTION ---------- //
	private SphLive sph;
	public static Particle[] particles;
	private static PrioQ pq;
	private Node root = new Node(0, numberOfParticles - 1, 0, 1, 0, 1);
	private JFrame jframe;
	private Dimension dim;
	private final int WIDTH = 1000, HEIGHT = 900, yFrame = 50;
	private final int displace = 50;
	private int width = 7, height = 7, pixel = 1;
	private List<Node> leafList = new ArrayList<Node>();
	private static int count = 0;
	public int numberOfLeafs = 0;
	private boolean isLeft = false;
	private int speed = 0;
	private int speedTreeSetUp = 100;
	private int bucketSize = 8;
	private boolean numbersTrue;
	public double gamma = 2.0;
	private double sigma = 0.15;
	private double alpha = 0.5, beta = 2 * alpha, eta = 0.007;
	static boolean firstSetUp = true;
	
	private double timeStep = 0.001;
	public static int numberOfParticles = 500;
	private static int numberOfSteps = 100000;
	
	
	public static void main(String[] args) {
		
		SphLive sph = new SphLive();
		sph.setJFrame(sph);
		sph.createParticles();
		sph.buildTreeDyn(particles, sph.root, 0, sph.bucketSize);
		sph.setLeafs(sph.root);
//		sph.printLeafList();
		sph.initialSearch();
		sph.setRho();
		sph.paintDyn("no");
		sph.run();	
}
	public int getNumberOfLeafs(){
		return this.numberOfLeafs; 
}
	
	private void run() {

		drift1(0);
		calcForces();
		
		for (int step = 0; step < numberOfSteps; step++) {
			drift1(timeStep / 2);
			calcForces();
			kick(timeStep);
			drift2(timeStep / 2);
				
				if (step % 10 == 0) {
					jframe.repaint();
					try {
						Thread.sleep(speed);
					} catch ( Exception e) {
						break;
					}
				}
			
			
//			if (step % 150 == 0) {
//				particles[numberOfParticles - 1] = new Particle(4., gaussRandom(sigma*2) + 0.5, gaussRandom(sigma*2) + 0.5);
//				particles[numberOfParticles - 1].v[0] = Math.random()*2. - 1;
//				particles[numberOfParticles - 1].v[1] = Math.random()*2. - 1;
//			}
		}
	}
	
	private void setJFrame(SphLive sph) {
		dim = Toolkit.getDefaultToolkit().getScreenSize();
		jframe = new JFrame("TreeDyn");
		jframe.setSize(WIDTH + 2*displace, HEIGHT + yFrame - 18 + 2*displace);
		jframe.setLocation(dim.width / 2 - (WIDTH + displace)/ 2, dim.height/ 2 - (HEIGHT+2*displace) / 2);
		jframe.setVisible(true);
//		jframe.setResizable(false);
		jframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		jframe.getContentPane().add(sph);
	}
	
	private static double square(double x) {
		return x*x;
	}
	
	private void createParticles() {
		particles = new Particle[numberOfParticles + 1]; //one Dummy-Particle (infinity coordinates)
		for (int i = 0; i < numberOfParticles; i++) {
			particles[i] = new Particle(1., gaussRandom(sigma) + 0.5, gaussRandom(sigma) + 0.5);
			
			particles[numberOfParticles] = new Particle(1, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
		}
	}
   
	private static double gaussRandom(double sigma) {
		double x1 = Math.random();
		double x2 = Math.random();
		
		double y = Math.sqrt(-2 * Math.log(x1)) * Math.cos(2 * Math.PI * x2) * sigma;
	
		while (y < -0.5 || y > 0.5) {
			double xNew1 = Math.random();
			double xNew2 = Math.random();				
			y = Math.sqrt(-2 * Math.log(xNew1)) * Math.cos(2 * Math.PI * xNew2) * sigma;
		}
		return y;
	}

	private static int partition(Particle[] particles, int i, int j, double s, int sel) {
		
		while (i < j){
			while ((i <= j) && (particles[i].r[sel] < s)) {
				i++;
			}
			while ((i <= j) && (particles[j].r[sel] >= s))  {
				j--;
			}
			if (i < j) {
				Particle tmp = particles[i];
				particles[i] = particles[j];
				particles[j] = tmp;
			}
		}
		return i;
	}
	

	private void buildTreeDyn(Particle[] particles, Node node, int sel, int bucketSize) {
		double s = (node.rmin[sel] + node.rmax[sel]) / 2;
				
		int k = partition(particles, node.start, node.end, s, sel);
		
		double[] rmin = new double[2];
		double[] rmax = new double[2];
		rmin[0] = node.rmin[0];
		rmax[0] = node.rmax[0];
		rmin[1] = node.rmin[1];
		rmax[1] = node.rmax[1];
		
		if ((k - 1) - node.start + 1 > bucketSize) {
			node.left = new Node(node.start, k - 1, rmin[0], rmax[0], rmin[1], rmax[1]);
			node.left.rmax[sel] = s;
			
			buildTreeDyn(particles, node.left, 1-sel, bucketSize);
		}
		
		else if ((k - 1) - node.start + 1 > 0) {
			node.left = new Node(node.start, k - 1, rmin[0], rmax[0], rmin[1], rmax[1]);
			node.left.rmax[sel] = s;
			
			for (int i = node.left.start; i <= node.left.end; i++) {
				node.left.balancePoint[0] += particles[i].r[0];
				node.left.balancePoint[1] += particles[i].r[1];
			}
			node.left.balancePoint[0] = node.left.balancePoint[0] / (node.left.nodeParticles);
			node.left.balancePoint[1] = node.left.balancePoint[1] / (node.left.nodeParticles);
		}
		
		if ((node.end - k + 1) > bucketSize) {
			node.right = new Node(k, node.end, rmin[0], rmax[0], rmin[1], rmax[1]);
			node.right.rmin[sel] = s;
			
			buildTreeDyn(particles, node.right, 1-sel, bucketSize);
		}
		else if (node.end - k + 1 > 0) {
			node.right = new Node(k, node.end, rmin[0], rmax[0], rmin[1], rmax[1]);
			node.right.rmin[sel] = s;
			
			for (int i = node.right.start; i <= node.right.end; i++) {
				node.right.balancePoint[0] += particles[i].r[0];
				node.right.balancePoint[1] += particles[i].r[1];
			}
			node.right.balancePoint[0] = node.right.balancePoint[0] / (node.right.nodeParticles);
			node.right.balancePoint[1] = node.right.balancePoint[1] / (node.right.nodeParticles);
		}

		if (node.right == null && node.left == null) {
			return;
		
		} else if (node.right == null) {
			node.balancePoint[0] = (node.left.nodeParticles * node.left.balancePoint[0]) / (node.left.nodeParticles);
			node.balancePoint[1] = (node.left.nodeParticles * node.left.balancePoint[1]) / (node.left.nodeParticles);
		} else if (node.left == null) {
			node.balancePoint[0] = (node.right.nodeParticles * node.right.balancePoint[0]) / (node.right.nodeParticles);
			node.balancePoint[1] = (node.right.nodeParticles * node.right.balancePoint[1]) / (node.right.nodeParticles);
		} else if (node.left != null && node.right != null) {
			node.balancePoint[0] = (node.left.nodeParticles * node.left.balancePoint[0] + node.right.nodeParticles * node.right.balancePoint[0]) / (node.left.nodeParticles + node.right.nodeParticles);
			node.balancePoint[1] = (node.left.nodeParticles * node.left.balancePoint[1] + node.right.nodeParticles * node.right.balancePoint[1]) / (node.left.nodeParticles + node.right.nodeParticles);
		}
	}
	
	private static void drift1(double timeStep) {
		drift2(timeStep);
		
		for (int i = 0; i < numberOfParticles; i++) {			
			particles[i].vPred[0] = particles[i].v[0] + particles[i].a[0] * timeStep;
			particles[i].vPred[1] = particles[i].v[1] + particles[i].a[1] * timeStep;
			
			particles[i].ePred = particles[i].e + particles[i].eDot * timeStep;
			if (particles[i].ePred < 0) {
				particles[i].ePred = 0;
			}
		}
	}
	
	private static void drift2(double timeStep) {
		for (int i = 0; i < numberOfParticles; i++) {
			particles[i].r[0] += particles[i].v[0] * timeStep;
			particles[i].r[1] += particles[i].v[1] * timeStep;
			
			for (int n = 0; n < 2; n++) {
				if (particles[i].r[n] < 0){
					particles[i].r[n] += 1;
				} else if (particles[i].r[n] >= 1) {
					particles[i].r[n] -= 1;
				}
			}
		}
	}
	
	private static void kick(double timeStep) {
		for (int i = 0; i < numberOfParticles; i++) {
			particles[i].v[0] += particles[i].a[0] * timeStep;
			particles[i].v[1] += particles[i].a[1] * timeStep;
			
			particles[i].e += particles[i].eDot * timeStep;
		}
	}
	
	private static double kernel(double r, double h, int D) {
		double sigma = 40 / (7 * Math.PI); // for D == 2
		double a = 0;
		if (0 <= r / h && r / h < 0.5) {
			a = 6 * Math.pow(r / h, 3) - 6 * square(r / h) + 1;
		} else if (0.5 <= r / h && r / h < 1) {
			a = 2 * Math.pow((1 - (r / h)), 3);
		} else {
			a = 0;
		}
		double omega = sigma / (Math.pow(h, D)) * a;
		return omega;
	}
	
	private void setRho () {
		double r = 0;
		double h = 0;
		for (int i = 0; i < numberOfParticles; i++) {
			particles[i].rho = 0;
			for (int k = 0; k < particles[i].pq.NN; k++) {
				int j = particles[i].pq.list[k];
				r = Math.sqrt(particles[i].pq.dist2[k]);
				h = particles[i].pq.getRadius();
				particles[i].rho += particles[j].mass * kernel(r, h, 2);
			}
		}
	}
	
	private void calcForces() {
		buildTreeDyn(particles, root, 0, bucketSize);
		for (int i = 0; i < numberOfParticles; i++) {
			
			particles[i].pq.init();
			search(particles[i].pq, root, i);
		}
		setRho();
		
		calcSound();
		NNSPHForces();
	}
	
	private void calcSound() {
		for (int i = 0; i < numberOfParticles; i++) {
			particles[i].c = Math.sqrt(gamma * (gamma - 1) * particles[i].ePred);
		}
	}
	
	private void NNSPHForces() {
		for (int i = 0; i < numberOfParticles; i++) {
			particles[i].P = square(particles[i].c) * particles[i].rho / gamma;
		}
		
		double r = 0;
		double h = 0;
		double pi = 0;
		double mu = 0;
		double[] vAB = new double[2];
		double[] rAB = new double[2];
		double rABDotvAB = 0;
		double delWdelr = 0;
		double gradW[] = new double[2];
		double sigma = 40 / (7 * Math.PI);
		double D = 2;
		double cAB = 0;
		double rhoAB = 0;
		
		for (int i = 0; i < numberOfParticles; i++) {
			particles[i].eDot = 0;
			particles[i].a[0] = 0;
			particles[i].a[1] = 0;
			
			for (int k = 0; k < particles[i].pq.NN; k++) {
				int j = particles[i].pq.list[k];
				if (i!=j) {
					
					r = Math.sqrt(particles[i].pq.dist2[k]);
					h = particles[i].pq.getRadius();
					
					//gradW
					if (r / h < 0.5) {
						delWdelr = 6 * sigma / Math.pow(h, D + 1) * ( 3 * square(r / h) - 2 * (r / h));
					} else if (r / h < 1) {
						delWdelr = 6 * sigma / Math.pow(h, D + 1) * (-1) * square(1 - r/h);
					} else {
						delWdelr = 0;
					}
					
					gradW[0] = (particles[i].pq.listdx[k]) / r * delWdelr;
					gradW[1] = (particles[i].pq.listdy[k]) / r * delWdelr;
					
					//eDot
					vAB[0] = particles[i].vPred[0] - particles[j].vPred[0]; vAB[1] = particles[i].vPred[1] - particles[j].vPred[1];
					rAB[0] = particles[i].pq.listdx[k]; rAB[1] = particles[i].pq.listdy[k];
					rABDotvAB = vAB[0] * rAB[0] + vAB[1] * rAB[1];
					mu =  h * rABDotvAB / (particles[i].pq.dist2[k] + eta*eta);
					cAB = 0.5 * (particles[i].c + particles[j].c);
					rhoAB = 0.5 * (particles[i].rho + particles[j].rho);
					if (rABDotvAB < 0){
						pi = ((-1) * alpha * cAB * mu + beta * square(mu)) / (rhoAB);
					} else {
						pi = 0;
					}
					
					//a[0] & a[1]
					for (int n = 0; n < 2; n++) {
						particles[i].a[n] += (-1) * particles[j].mass * (particles[i].P / square(particles[i].rho) + particles[j].P / square(particles[j].rho) + pi) * gradW[n];
					}
					
					particles[i].eDot += particles[j].mass * (vAB[0] * gradW[0] + vAB[1] * gradW[1]);
				}
			}
			particles[i].eDot *= particles[i].P / square(particles[i].rho);
			
		}
	}
		
		
	private static void search(PrioQ pq, Node c, int i) {
		searchPeriodic(pq, c, i, new double[]{0,0});
		searchPeriodic(pq, c, i, new double[]{-1,-1});
		searchPeriodic(pq, c, i, new double[]{0,-1});
		searchPeriodic(pq, c, i, new double[]{1,-1});
		searchPeriodic(pq, c, i, new double[]{1,0});
		searchPeriodic(pq, c, i, new double[]{1,1});
		searchPeriodic(pq, c, i, new double[]{0,1});
		searchPeriodic(pq, c, i, new double[]{-1,1});
		searchPeriodic(pq, c, i, new double[]{-1,0});
	}
	
	private static void searchPeriodic(PrioQ pq, Node c, int i, double[] offset) {
		double d2left;
		double d2right;
		double d2;
		if (c == null) {
			return;
		} if (c.right != null && c.left != null) {
			d2left = square((c.left.balancePoint[0] + offset[0] - particles[i].r[0])) + square((c.left.balancePoint[1] + offset[1] - particles[i].r[1]));
			d2right = square((c.right.balancePoint[0] + offset[0] - particles[i].r[0])) + square((c.right.balancePoint[1] + offset[1] - particles[i].r[1]));
			
			if (d2left < d2right) {
				if (pq.mindist2(c.left, particles[i], offset) < (pq.dist2[pq.index])) {
					searchPeriodic(pq, c.left, i, offset);
				}
				if (pq.mindist2(c.right, particles[i], offset) < (pq.dist2[pq.index])){
					searchPeriodic(pq, c.right, i, offset);
				}
			}else { // d2right < d2left
				if (pq.mindist2(c.right, particles[i], offset) < (pq.dist2[pq.index])) {
					searchPeriodic(pq, c.right, i, offset);
				}
				if (pq.mindist2(c.left, particles[i], offset) < (pq.dist2[pq.index])){
					searchPeriodic(pq, c.left, i, offset);
				}
			}
		} else if (c.left != null) {
			searchPeriodic (pq, c.left, i, offset);
		} else if (c.right != null) {
			searchPeriodic (pq, c.right, i, offset);
		} else { //leaf cell
			for (int j = c.start; j <= c.end; j++) {
				double dx = (particles[i].r[0] - (particles[j].r[0] + offset[0]));
				double dy = (particles[i].r[1] - (particles[j].r[1] + offset[1]));
				pq.replace(dx, dy, j);
			}
		}
	}
	
	private void initialSearch() {
		for (int i = 0; i < particles.length; i++) {
			pq = new PrioQ();
			pq.init();
			search(pq, root, i);
			particles[i].pq = pq;
		}
	}
	
	private void setLeafs(Node node) {
		if (node.left == null && node.right == null) {
			leafList.add(node);
			numberOfLeafs++;
			leafList.get(numberOfLeafs - 1).setLeftTrue(isLeft);
		}
		
		if (node.left != null){
			isLeft = true;
			setLeafs(node.left);
		}
		
		if (node.right != null){
			isLeft = false;
			setLeafs(node.right);
		}
	}
	
	private void paintDyn(String String2) {
		firstSetUp = true;
		if (String2.equals("yes")) {
			numbersTrue = true;
		} else if (String2.equals("no")) {
			numbersTrue = false;
		} else {
			numbersTrue = false;
		}
		
		while (count < numberOfLeafs) {
			
			jframe.repaint(); //calls paintComponent();
			
			count++;
			try {
				if (count==numberOfLeafs){
					speedTreeSetUp = 3000;
				}
				Thread.sleep(speedTreeSetUp);
			} catch ( Exception e) {
				break;
			}
		}
		count=0;
		firstSetUp=false;
		
	}
	
	protected void paintComponent(Graphics g) {
		double rhoFactor = 3;
		double rhomax = numberOfParticles * rhoFactor;
		if (firstSetUp==false){
			//Afterwards, painting SPH particles simulation
			int xs_max = WIDTH, xs;
			int ys_max = HEIGHT, ys;
			width=7;
			height=7;
			double scale = 70;
			double dx, dy;
			
			
			g.setColor(Color.white);
			g.fillRect(displace, displace, WIDTH, HEIGHT);
	//		System.out.println("Leaf number: " + sph.getNumberOfLeafs());
			for (int i = 0; i < numberOfParticles; i++) {
				xs = (int)(particles[i].r[0] * xs_max) + displace;
				ys = (int)(particles[i].r[1] * ys_max) + displace;
				double cv = particles[i].rho/rhomax;
				if (cv < 0) {
					cv = 0;
				} else if (cv > 1) {
					cv = 1;
				}
				
				g.setColor(Rainbow(1-cv));
				
				dx = particles[i].v[0] * timeStep * scale;
				dy = particles[i].v[1] * timeStep * scale;
				
				g.fillOval(xs - width / 2, ys - height / 2, width, height);
				g.drawLine(xs, ys, (int)(xs_max * (particles[i].r[0] + dx)) + displace, (int)(ys_max * (particles[i].r[1] + dy ))+ displace);
				int headMinus = 4;
				g.fillOval((int)(xs_max * (particles[i].r[0] + dx)-(width-headMinus)/2) + displace, (int)(ys_max * (particles[i].r[1] + dy)-(width-headMinus)/2) + displace, width-headMinus, height-headMinus);
			}
		} else if (firstSetUp == true){
			//First build up, painting tree structure
			super.paintComponent(g);
			int xs_max = WIDTH, xs;
			int ys_max = HEIGHT, ys;
			
			g.setColor(Color.white);
			g.fillRect(displace, displace, WIDTH, HEIGHT);
			for (int i = 0; i < count; i++) {
				
				if (leafList.get(i).getLeftTrue() == true) {
					g.setColor(new Color((float)0, (float)1, (float)0, (float)0.3));
				} else {
					g.setColor(Color.red);
					g.setColor(new Color((float)1, (float)0, (float)0, (float)0.3));
				}
				
				g.fillRect((int)(leafList.get(i).rmin[0] * xs_max) + pixel + displace, (int)(leafList.get(i).rmin[1]* ys_max) + pixel+ displace, (int)((leafList.get(i).rmax[0] - leafList.get(i).rmin[0]) * xs_max) - 2*pixel, (int)((leafList.get(i).rmax[1] - leafList.get(i).rmin[1]) * ys_max) - 2*pixel);
				
				if (numbersTrue) {
					g.setColor(Color.black);
					g.drawString(String.valueOf(i), (int)(leafList.get(i).rmin[0] * xs_max) + 2, (int)(leafList.get(i).rmin[1]* ys_max) + 13);
				}
				
//				System.out.println("Balance-Point of leaf cell " + i + ": (x,y) (" + leafList.get(i).balancePoint[0] + ", " + leafList.get(i).balancePoint[1] + ")" );
				
				//draw particles inside leaf cells
				for (Particle particle : particles) {
					if (particle.r[0] > leafList.get(i).rmin[0] && particle.r[0] < leafList.get(i).rmax[0] 
							&& particle.r[1] > leafList.get(i).rmin[1] && particle.r[1] < leafList.get(i).rmax[1]) {
						xs = (int)(particle.r[0] * xs_max)+ displace;
						ys = (int)(particle.r[1] * ys_max)+ displace;
						double cv = particle.rho/rhomax;
						if (cv < 0) {
							cv = 0;
						} else if (cv > 1) {
							cv = 1;
						}
						
						g.setColor(Rainbow(1-cv));
						
//						float alpha = (float) 1.0;
//						float r,green,b;
//						r = 0; green = 0; b = 1;
//						g.setColor(new Color(r,green,b,alpha));
						g.fillOval(xs - width / 2, ys - height / 2, width, height);
					}
				}
				
				
				
//				//draw circle for 32 nearest neighbors
//				for (Particle particle : particles) {
//					xs = (int)(particle.r[0] * xs_max);
//					ys = (int)(particle.r[1] * ys_max);
//					g.setColor(Color.black);
//					g.drawOval(xs - (int)(particle.pq.getRadius() * xs_max), ys - (int)(particle.pq.getRadius() * ys_max), (int)(2 * particle.pq.getRadius() * xs_max), (int)(2 * particle.pq.getRadius() * ys_max));
//				System.out.println("ARIIVED END");
//				}
			}
			
			
		}
		
		
	}
	
	
	//Red-Green Leafs
	public void paintComponent2(Graphics g) {
		super.paintComponent(g);
		int xs_max = WIDTH, xs;
		int ys_max = HEIGHT, ys;
		
		g.setColor(Color.white);
		g.fillRect(0, 0, WIDTH, HEIGHT);
		for (int i = 0; i < count; i++) {
			if (leafList.get(i).getLeftTrue() == true) {
				g.setColor(new Color((float)0, (float)1, (float)0, (float)0.7));
			} else {
				g.setColor(Color.red);
				g.setColor(new Color((float)1, (float)0, (float)0, (float)0.7));
			}
			
			g.fillRect((int)(leafList.get(i).rmin[0] * xs_max) + pixel + displace, (int)(leafList.get(i).rmin[1]* ys_max) + pixel+ displace, (int)((leafList.get(i).rmax[0] - leafList.get(i).rmin[0]) * xs_max) - 2*pixel, (int)((leafList.get(i).rmax[1] - leafList.get(i).rmin[1]) * ys_max) - 2*pixel);
			
			if (numbersTrue) {
				g.setColor(Color.black);
				g.drawString(String.valueOf(i), (int)(leafList.get(i).rmin[0] * xs_max) + 2, (int)(leafList.get(i).rmin[1]* ys_max) + 13);
			}
			
//			System.out.println("Balance-Point of leaf cell " + i + ": (x,y) (" + leafList.get(i).balancePoint[0] + ", " + leafList.get(i).balancePoint[1] + ")" );
			
			//draw particles inside leaf cells
			for (Particle particle : particles) {
				if (particle.r[0] > leafList.get(i).rmin[0] && particle.r[0] < leafList.get(i).rmax[0] 
						&& particle.r[1] > leafList.get(i).rmin[1] && particle.r[1] < leafList.get(i).rmax[1]) {
					xs = (int)(particle.r[0] * xs_max)+ displace;
					ys = (int)(particle.r[1] * ys_max)+ displace;
					float alpha = (float) 1.0;
					float r,green,b;
					r = 0; green = 0; b = 1;
					g.setColor(new Color(r,green,b,alpha));
					g.drawOval(xs - width / 2, ys - height / 2, width, height);
				}
			}
			
			
//			//draw circle for 32 nearest neighbors
//			for (Particle particle : particles) {
//				xs = (int)(particle.r[0] * xs_max);
//				ys = (int)(particle.r[1] * ys_max);
//				g.setColor(Color.black);
//				g.drawOval(xs - (int)(particle.pq.getRadius() * xs_max), ys - (int)(particle.pq.getRadius() * ys_max), (int)(2 * particle.pq.getRadius() * xs_max), (int)(2 * particle.pq.getRadius() * ys_max));
//			System.out.println("ARIIVED END");
//			}
		}
	}
	
//	public void paintComponent(Graphics g) {
//	
//	paintTree(root, g, WIDTH / 2, 0, 200, 20);
//}
//
//private void paintTree(Node node, Graphics g, int x, int y, int plusX, int plusY) {
//	g.fillOval(x, y, width, height);
//	if (node.left == null && node.right == null) {
//		
//	}
//	
//	if (node.left != null){
//		g.setColor(Color.green);
//		g.drawLine(x, y, x - plusX, y + plusY);
//		paintTree(node.left, g, x - plusX, y + plusY, plusX / 2, plusY);
//	}
//	
//	if (node.right != null){
//		g.setColor(Color.red);
//		g.drawLine(x, y, x + plusX, y + plusY);
//		paintTree(node.right, g, x + plusX, y + plusY, plusX / 2, plusY);
//	}
//}
	private static Color Rainbow(double x) {
		int r,g,b;
		
		if ( 0 <= x && x < 0.2) {
			r = 255;
			b = 0;
			g = (int)(x / 0.2 * 256);
		} else if (0.2 <= x && x < 0.4) {
			r = (int)((1-(x - 0.2) / 0.2) * 255);
			b = 0;
			g = 255;
		} else if (0.4 <= x && x < 0.6) {
			r = 0;
			b = (int)((x - 0.4) / 0.2 * 255);
			g = 255;
		} else if (0.6 <= x && x < 0.8) {
			r = 0;
			b = 255;
			g = (int)((1-(x - 0.6) / 0.2) * 255);
		} else if (0.8 <= x && x <= 1.0) {
			r = (int)((x - 0.8) / 0.2 * 255);
			b = 255;
			g = 0;
		} else if (x > 1) {
			r = g = b = 255; 
			
		} else {
			r = 255;
			g = 0;
			b = 255;
		}
		return new Color(r,g,b);
	}
	
	@SuppressWarnings("unused")
	private static void printParticleList() {
		System.out.println("Particle-List: " + particles.length + " elements.");
		for (int i = 0; i < numberOfParticles; i++) {
			System.out.println("-----------------------------");
			System.out.println("Index of Particle, number: " + i);
			System.out.println("x-value: " + particles[i].r[0]);
			System.out.println("y-value: " + particles[i].r[1]);
			System.out.println("-----------------------------");
			System.out.println();
		}
	}
	
	@SuppressWarnings("unused")
	private void printLeafList() {
		System.out.println("Leaf-List: " + leafList.size() + " elements.");
		System.out.println("Number of Leafs: " + numberOfLeafs);
		for (int i = 0; i < leafList.size(); i++) {
			System.out.println("-----------------------------");
			System.out.println("Index of Leaf: " + i);
			if (leafList.get(i).getLeftTrue()){
				System.out.println("LEFT");
			} else {
				System.out.println("RIGHT");
			}
			System.out.println("x-range: " + leafList.get(i).rmin[0] + " to " + leafList.get(i).rmax[0]);
			System.out.println("y-range: " + leafList.get(i).rmin[1] + " to " + leafList.get(i).rmax[1]);
			System.out.println("-----------------------------");
			System.out.println();
		}
	}
	
	private void printNode(Node node, int i) {
		System.out.println("Level: " + i);
		System.out.println("Balance Points (x, y): (" + node.balancePoint[0] + ", " + node.balancePoint[1] + ")");
		if (node.left != null) {
			printNode(node.left, i + 1);
		}
		if (node.right != null) {
			printNode(node.right, i + 1);
		}
	}
	
	private static void printRho() {
		for (Particle particle : particles) {
			System.out.println("Rho: " + particle.rho);
		}
	}
		
}





