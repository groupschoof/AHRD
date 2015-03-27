package ahrd.controller;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

public class AddEuclideanDistsToPathLogTable {

	public static final String PATH_LOG_TABLE_SEP_CHAR = "\t";
	public static final String DB_WEIGHT_COLS_ARG_SEP_CHAR = ",";
	public static final double NORMALIZE_DB_WEIGHT_DENOMINATOR = 100;
	public static final String PATH_LOG_EUCLIDEAN_DIST_COL_NAME = "Euclidean.Distance.2.Curr.Acptd";
	public static final int REJECTED_WORSE_PERFORMING_PARAMS = 0;

	public static void main(String[] args) throws IOException {
		System.out
				.println("Usage: java -Xmx2g -cp ahrd.jar ahrd.controller.AddEuclideanDistsToPathLogTable path_log_table.tsv startColumn stopColumn acceptedColumn temperatureColumn dbWeightsColInds(e.g. 5,8,11) extended_path_log_table_out.tsv");
		BufferedReader tblIn = null;
		BufferedWriter tblOut = null;
		int startCol = Integer.parseInt(args[1]);
		// Arrays.copyOfRange excludes the 'until' column, hence add 1:
		int stopCol = Integer.parseInt(args[2]) + 1;
		int acceptedCol = Integer.parseInt(args[3]);
		int tempeCol = Integer.parseInt(args[4]);
		int[] dbWghtsCols = parseDatabaseWeightColumnArg(args[5], startCol);
		try {
			// In- and Output:
			tblIn = new BufferedReader(new FileReader(args[0]));
			tblOut = new BufferedWriter(new FileWriter(args[6]));
			// Iterated lines:
			String header, line;
			long currTemp, lastTemp;
			if ((header = tblIn.readLine()) != null) {
				if ((line = tblIn.readLine()) != null) {
					// Write first lines into output:
					writeHeaderToOutputPathLogTable(tblOut, header, line);
					// Initialize first parameter set in path log table:
					double[] currAcpt = parseCurrentPathLogLine(line, startCol,
							stopCol);
					lastTemp = parseTemperature(line, tempeCol);
					double[] currEval;
					while ((line = tblIn.readLine()) != null) {
						currEval = parseCurrentPathLogLine(line, startCol,
								stopCol);
						currTemp = parseTemperature(line, tempeCol);
						double euclDist = currTemp <= lastTemp ? measureEuclideanDist(
								normalizeDatabaseWeights(currAcpt, dbWghtsCols),
								normalizeDatabaseWeights(currEval, dbWghtsCols))
								: 0;
						appendLineToOutputPathLogTable(tblOut, line, euclDist);
						// Set currently accepted parameters, if the currently
						// evaluated ones have been accepted during simulated
						// annealing:
						if (isAccepted(line, acceptedCol))
							currAcpt = currEval;
						// Check using temperature, if not a new batch of path
						// logs has started:
						lastTemp = currTemp;
					}
				}
			}
		} finally {
			tblIn.close();
			tblOut.close();
		}
	}

	public static boolean isAccepted(String line, int acceptedCol) {
		int acpt = Integer
				.parseInt(line.split(PATH_LOG_TABLE_SEP_CHAR)[acceptedCol]);
		return acpt != REJECTED_WORSE_PERFORMING_PARAMS;
	}

	public static void writeHeaderToOutputPathLogTable(BufferedWriter out,
			String header, String firstParams) throws IOException {
		out.write(header + PATH_LOG_TABLE_SEP_CHAR
				+ PATH_LOG_EUCLIDEAN_DIST_COL_NAME);
		out.newLine();
		out.write(firstParams + PATH_LOG_TABLE_SEP_CHAR + "0.0");
		out.newLine();
	}

	public static void appendLineToOutputPathLogTable(BufferedWriter out,
			String srcLine, double euclDist) throws IOException {
		out.write(srcLine + PATH_LOG_TABLE_SEP_CHAR + euclDist);
		out.newLine();
	}

	public static double[] parseCurrentPathLogLine(String line, int startCol,
			int stopCol) {
		double[] params = new double[stopCol - startCol];
		int i = 0;
		for (String iParam : Arrays.copyOfRange(
				line.split(PATH_LOG_TABLE_SEP_CHAR), startCol, stopCol)) {
			params[i] = Double.parseDouble(iParam);
			i++;
		}
		return params;
	}

	public static int[] parseDatabaseWeightColumnArg(String dbWghtsArg,
			int startCol) {
		String[] dbWeightColsArgs = dbWghtsArg
				.split(DB_WEIGHT_COLS_ARG_SEP_CHAR);
		int[] dbWghtsCols = new int[dbWeightColsArgs.length];
		for (int i = 0; i < dbWeightColsArgs.length; i++) {
			dbWghtsCols[i] = Integer.parseInt(dbWeightColsArgs[i]) - startCol;
		}
		return dbWghtsCols;
	}

	public static double[] normalizeDatabaseWeights(double[] params,
			int[] dbWghtsCols) {
		double[] normPars = Arrays.copyOf(params, params.length);
		for (int i : dbWghtsCols) {
			normPars[i] = params[i] / NORMALIZE_DB_WEIGHT_DENOMINATOR;
		}
		return normPars;
	}

	public static double measureEuclideanDist(double[] x, double[] y) {
		if (x.length != y.length) {
			throw new IllegalArgumentException(
					"Cannot measure euclidean distance between two double[] of unequal length.");
		}
		double dist = 0;
		for (int i = 0; i < x.length; i++) {
			dist += Math.pow(x[i] - y[i], 2);
		}
		return Math.sqrt(dist);
	}

	public static long parseTemperature(String pathLogLine,
			int temparatureColumnIndex) {
		return Integer
				.parseInt(pathLogLine.split(PATH_LOG_TABLE_SEP_CHAR)[temparatureColumnIndex]);
	}
}
