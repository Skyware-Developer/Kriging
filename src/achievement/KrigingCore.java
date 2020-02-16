package achievement;

import Jama.Matrix;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.stream.Collectors;
import java.text.DecimalFormat;

import bean.DEM_Point;
import in_and_out.Input;
import achievement.meshing;

public class KrigingCore {

	/**
	 * Kriging
	 * 
	 * @return
	 */
	// Semi-variance function, also called variogram. Spherical model formula
	// (spherical model is usually used when interpolating meteorological elements)
	// Initial estimated range a = 130.000 Abutment value is 850.000 Nugget constant
	// is 500.000
	public static double function(double h) {
		// Initialize variogram value
		double r = 0;
		// h is the step size
		if (h == 0) {
			r = 500.000;
		} else if (h > 0 && h <= 130) {
			r = 500.000 + 13.500 * (((3 * h) / (2 * 130)) - ((h * h * h) / (2 * 130 * 130 * 130)));
			// r = 0.000 + 13.500 * ((3 * h * h * h * h) / (4 * 30 * 30 * 30 * 30));
		} else {
			// r = C0+C Abutment value
			r = 850.000;
		}
		return r;

	}
	// Variation function using Gaussian model
	/*
	 * public double function(double h) { double r = 0; if(h == 0) r = 0; else if(h
	 * > 0 && h <= 11 ) { r = 0.0 + 0.06 * (1.0 - Math.exp(-h/11) ); } else r =
	 * 3.202;
	 * 
	 * return r; }
	 */

	public static ArrayList<DEM_Point> gridData()
	// public static void gridData()
	{
		ArrayList<DEM_Point> demPoints = meshing.border();// Grid point array
		ArrayList<ArrayList<DEM_Point>> sum = Input.readFile();// Discrete point data
		// System.out.println(demPoints);

		int count = 0;
		int XLength = sum.size();
		// System.out.println( XLength ) ;
		ArrayList<DEM_Point> dem = new ArrayList<DEM_Point>();// Grid point data after interpolation
		double[] Xs = new double[XLength];// Array of discrete points X
		double[] Ys = new double[XLength];// Array of discrete points Y
		double[] Zs = new double[XLength];// Array of discrete points Z
		double[] rData = new double[XLength + 1];// Distance between grid points and discrete points
		double[][] ZD = new double[XLength + 1][XLength + 1];// Difference of discrete points Z
		double[][] data = new double[XLength + 1][XLength + 1];// Distance between discrete points
		List<DEM_Point> points = null;

		// System.out.println(sum.size());
		// Take discrete points and divide them into an array Xs[],Ys[],Zs[]
		for (int i = 0; i < sum.size(); i++) {
			points = new ArrayList<DEM_Point>(sum.get(i));
			// System.out.println(i+" "+sum.get(i));
			List<Double> X = points.stream().map(DEM_Point::getX).collect(Collectors.toList());
			// System.out.println(X);
			Xs[i] = X.get(0);
			// System.out.println(Xs[i]);
			List<Double> Y = points.stream().map(DEM_Point::getY).collect(Collectors.toList());
			Ys[i] = Y.get(0);
			List<Double> Z = points.stream().map(DEM_Point::getZ).collect(Collectors.toList());
			Zs[i] = Z.get(0);
		}

		/*
		 * for(int i = 0; i < Xs.length; i++) { System.out.println(Xs[i]); }
		 */

		// Find the variation function, which is the coefficient matrix
		for (int i = 0; i <= XLength; i++) {
			for (int j = 0; j <= XLength; j++) {
				if ((i < XLength) && (j < XLength)) {
					if (i != j) {
						double TempDis = 0;

						TempDis = Math.sqrt(((Xs[i] - Xs[j]) * (Xs[i] - Xs[j])) + ((Ys[i] - Ys[j]) * (Ys[i] - Ys[j])));
						ZD[i][j] = Math.abs(Zs[i] - Zs[j]);// Find the absolute value

						// System.out.println(TempDis);
						data[i][j] = function(TempDis);
						/*
						 * data[i][j] = Math.sqrt( ( (Xs[i]-Xs[j]) * (Xs[i]-Xs[j]) ) + ( (Ys[i]-Ys[j]) *
						 * (Ys[i]-Ys[j]) ) );
						 */
					} else {
						data[i][j] = 0;
					}
				} else if ((i == XLength) || (j == XLength)) {
					data[i][j] = 1;
				} else if ((i == XLength) && (j == XLength)) {
					data[i][j] = 0;
				}
			}
		}

		/*
		 * for(int i = 0; i <XLength ; i++ ) { for(int j = 0; j < XLength ; j++ ) {
		 * System.out.println(data[i][j] + "  " + i + "  " + j); } }
		 */

		Matrix m = new Matrix(data);// Get matrix m from two-dimensional array data
		Matrix M = m.inverse(); // Matrix M obtained by inverting matrix m

		while (count < demPoints.size()) {
			DEM_Point p = demPoints.get(count++);
			double xp = p.getX();// Grid Point X
			double yp = p.getY();// Grid Point Y
			double zp = 0;// 格网点Z
			// System.out.println(p);
			for (int k = 0; k <= XLength; k++) {
				if (k < XLength) {
					double TempDis2 = 0;

					TempDis2 = Math.sqrt((xp - Xs[k]) * (xp - Xs[k]) + (yp - Ys[k]) * (yp - Ys[k]));
					rData[k] = function(TempDis2);
				} else
					rData[k] = 1;
				// System.out.println(rData[k]);
			}
			// System.out.println(" ------------------------------------------");
			double[] rData2 = new double[XLength + 1];// Weighting factor λi
			for (int k = 0; k < (XLength + 1); k++) {
				for (int k2 = 0; k2 < (XLength + 1); k2++) {
					rData2[k] += M.get(k, k2) * rData[k2];
				}
				// System.out.println(rData2[k]);
			}

			// System.out.println(" ------------------------------------------");
			for (int k = 0; k < XLength; k++)// Get estimated value
			{
				zp += Zs[k] * rData2[k];
				// System.out.println(rData2[k]);
				// System.out.println(Zs[k] +" "+ rData2[k]+" "+ count);
			}

			DecimalFormat df = new DecimalFormat("###.000");
			double Zp = Double.parseDouble(df.format(zp));// Three decimal places
			// System.out.println(zp+" "+ count);
			dem.add(new DEM_Point(xp, yp, Zp));
		}
		// List<Double> X =
		// points.stream().map(DEM_Point::getX).collect(Collectors.toList());
		System.out.println(dem.stream().map(DEM_Point::getX));
		return dem;
	}

	public static void main(String[] args) { // TODO Auto-generated method stub
		System.out.println("开始执行操作");
		System.out.println("Start execution");
		ArrayList<DEM_Point> array = gridData();
		for (Iterator<DEM_Point> iterator = array.iterator(); iterator.hasNext();) {
			DEM_Point dem_Point = (DEM_Point) iterator.next();
			System.out.println(dem_Point.toString());
		}
		System.out.println(array.size());
		System.out.println("完成操作");
		System.out.println("Complete operation");

	}

}
