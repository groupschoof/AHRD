package ahrd.controller;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

public class AddEuclidieanDistsToPathLogTable {

	public static final String PATH_LOG_TABLE_SEP_CHAR = "\t";
	public static final String DB_WEIGHT_COLS_ARG_SEP_CHAR = ",";
	public static final double NORMALIZE_DB_WEIGHT_DENOMINATOR = 100;
	public static final String PATH_LOG_EUCLIDEAN_DIST_COL_NAME = "Euclidean.Distance.2.Curr.Acptd";

	public static void main(String[] args) throws IOException {
		System.out
				.println("Usage: java -Xmx2g -cp ahrd.jar ahrd.controller.AddEuclideanDistsToPathLogTable path_log_table.tsv startColumn stopColumn dbWeightsColInds(e.g. 5,8,11) extended_path_log_table_out.tsv");
		BufferedReader tblIn = null;
		BufferedWriter tblOut = null;
		int startCol = Integer.parseInt(args[1]);
		int stopCol = Integer.parseInt(args[2]);
		int[] dbWghtsCols = parseDatabaseWeightColumnArg(args[3]);
		try {
			// In- and Output:
			tblIn = new BufferedReader(new FileReader(args[0]));
			tblOut = new BufferedWriter(new FileWriter(args[4]));
			// Iterated lines:
			String header, line;
			if ((header = tblIn.readLine()) != null) {
				if ((line = tblIn.readLine()) != null) {
					// Write first lines into output:
					writeHeaderToOutputPathLogTable(tblOut, header, line);
					// Initialize first parameter set in path log table:
					double[] currAcpt = normalizeDatabaseWeights(
							parseCurrentPathLogLine(line, startCol, stopCol),
							dbWghtsCols);
					double[] currEval;
					while ((line = tblIn.readLine()) != null) {
						currEval = normalizeDatabaseWeights(
								parseCurrentPathLogLine(line, startCol, stopCol),
								dbWghtsCols);
						appendLineToOutputPathLogTable(tblOut, line, currEval,
								measureEuclideanDist(currAcpt, currEval),
								startCol);
					}
				}
			}
		} finally {
			tblIn.close();
			tblOut.close();
		}
	}

	private static void writeHeaderToOutputPathLogTable(BufferedWriter out,
			String header, String firstParams) throws IOException {
		out.write(header + PATH_LOG_TABLE_SEP_CHAR
				+ PATH_LOG_EUCLIDEAN_DIST_COL_NAME);
		out.newLine();
		out.write(firstParams + PATH_LOG_TABLE_SEP_CHAR + "0.0");
		out.newLine();
	}

	private static void appendLineToOutputPathLogTable(BufferedWriter out,
			String srcLine, double[] currEvalParams, double euclDist,
			int startCol) throws IOException {
		String[] cols = Arrays.copyOfRange(
				srcLine.split(PATH_LOG_TABLE_SEP_CHAR), 0, startCol - 1);
		StringBuilder builder = new StringBuilder();
		for (String iCol : cols) {
			builder.append(iCol);
			builder.append(PATH_LOG_TABLE_SEP_CHAR);
		}
		for (int i = 0; i < currEvalParams.length; i++) {
			builder.append(currEvalParams[i]);
			builder.append(PATH_LOG_TABLE_SEP_CHAR);
		}
		builder.append(euclDist);
		out.write(builder.toString());
		out.newLine();
	}

	private static double[] parseCurrentPathLogLine(String line, int startCol,
			int stopCol) {
		double[] params = new double[stopCol - startCol];
		int i = 0;
		for (String iParam : Arrays.copyOfRange(
				line.split(PATH_LOG_TABLE_SEP_CHAR), startCol, stopCol)) {
			params[i] = Integer.parseInt(iParam);
			i++;
		}
		return params;
	}

	private static double[] parseDatabaseWeightColumnArg(String dbWghtsArg) {
		String[] dbWeightColsArgs = dbWghtsArg
				.split(DB_WEIGHT_COLS_ARG_SEP_CHAR);
		double[] dbWghtsCols = new double[dbWeightColsArgs.length];
		for (int i = 0; i < dbWeightColsArgs.length; i++) {
			dbWghtsCols[i] = Double.parseDouble(dbWeightColsArgs[i]);
		}
		return dbWghtsCols;
	}

	private static double[] normalizeDatabaseWeights(double[] params,
			int[] dbWghtsCols) {
		double[] normPars = Arrays.copyOf(params, params.length);
		for (int i : dbWghtsCols) {
			normPars[i] = params[i] / NORMALIZE_DB_WEIGHT_DENOMINATOR;
		}
		return normPars;
	}

	private static double measureEuclideanDist(double[] x, double[] y) {
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
}
