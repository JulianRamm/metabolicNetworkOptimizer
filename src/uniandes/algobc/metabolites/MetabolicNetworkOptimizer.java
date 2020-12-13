package uniandes.algobc.metabolites;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class MetabolicNetworkOptimizer {

	private MetabolicNetwork inputNetwork;
	private Set<Metabolite> protectedMetabolites ;
	private Set<Reaction> protectedReactions;
	private Map<Inequality, Boolean> protectedFunctionsAndPhenotypes;
	private int minimunDegreesOfFreedom;
	private int minimumNumberOfReactions;

	public MetabolicNetworkOptimizer(MetabolicNetwork inputNetwork, Set<Metabolite> protectedMetabolites,
			Set<Reaction> protectedReactions, Map<Inequality, Boolean> protectedFunctionsAndPhenotypes,
			int minimunDegreesOfFreedom, int minimumNumberOfReactions) {
		this.inputNetwork = inputNetwork;
		this.protectedMetabolites = protectedMetabolites;
		this.protectedReactions = protectedReactions;
		this.protectedFunctionsAndPhenotypes = protectedFunctionsAndPhenotypes;
		this.minimunDegreesOfFreedom = minimunDegreesOfFreedom;
		this.minimumNumberOfReactions = minimumNumberOfReactions;
	}

	public MetabolicNetwork optimize() {
		if (preprocessing()) {
			List<Reaction> removables = removables();
			int dof = calculateDegreesOfFreedom();
			System.out.println(dof);
			while (dof > minimunDegreesOfFreedom && notEmpty(removables)
					&& inputNetwork.getReactionsList().size() > minimumNumberOfReactions) {
				fva(removables);
				boolean success = false;
				while (success == false && notEmpty(removables)) {
					Reaction cand = candidateReactForRemoval(removables);
					removeReaction(cand, removables);
					success = checkProtectedFunctions();

					if (success == false) {
						reinsert(cand);
					}
				}
			}
		}
		return inputNetwork;
	}

	public boolean preprocessing() {
		return checkProtectedFunctions();
	}

	public boolean checkProtectedFunctions() {
		List<Reaction> reactions = inputNetwork.getReactionsList();
		boolean res = true;
		for (Map.Entry<Inequality, Boolean> e : protectedFunctionsAndPhenotypes.entrySet()) {
			boolean val = true;
			for (int i = 0; i < reactions.size(); i++) {
				Inequality in = new Inequality();
				Reaction r = reactions.get(i);
				Inequality actual = e.getKey();
				if (actual.getRHSValues().length > 0) {
					in.createInequality(actual.getLHSConstant(), actual.getRHSConstant(), 0.3, r.getUpperBoundFlux(),
							r.getUpperBoundFlux(), r.getUpperBoundFlux(), r.getUpperBoundFlux());
				}
				in.createInequality(actual.getLHSConstant(), actual.getRHSConstant(), 0.3, r.getUpperBoundFlux(),
						r.getUpperBoundFlux());
				val = val ||( actual.checkFeasibility() && in.checkFeasibility());
			}
			e.setValue(val);
		}
		for (Map.Entry<Inequality, Boolean> e : protectedFunctionsAndPhenotypes.entrySet()) {
			res = res && e.getValue();
		}
		return res;
	}

	public void fva(List<Reaction> removables) {
		for (int i = 0; i < removables.size(); i++) {
			if (checkRemovability(removables.get(i)) == false) {
				removables.remove(i);
			}
		}
	}

	public boolean checkRemovability(Reaction r) {
		boolean res = false;
		double[] act = { 0, 0 };
		for (Map.Entry<Inequality, Boolean> e : protectedFunctionsAndPhenotypes.entrySet()) {
			act = e.getKey().getFeasibilityRange(r.getLowerBoundFlux(), r.getUpperBoundFlux(), 0.2);
			res = res || ((act[0] <= 0 && act[1] <= 0) || (act[0] >= 0 && act[1] >= 0));
		}
		return res;
	}

	public double[] calculateFluxRange(Reaction r) {
		double res[] = { 0, 0 };
		ArrayList<double[]> ranges = new ArrayList<>();
		for (Map.Entry<Inequality, Boolean> e : protectedFunctionsAndPhenotypes.entrySet()) {
			ranges.add(e.getKey().getFeasibilityRange(r.getLowerBoundFlux(), r.getUpperBoundFlux(), 0.2));
		}
		for (int i = 0; i < ranges.size(); i++) {
			double[] act = ranges.get(i);
			for (int j = 0; j < ranges.size(); j++) {
				double[] compare = ranges.get(j);
				if (compare[0] >= act[0] || compare[1] <= act[1]) {
					act[0] = Math.min(compare[0], act[0]);
					act[1] = Math.max(compare[1], act[1]);
				}
			}
			if (act[0] >= res[0] || act[1] <= res[1]) {
				res[0] = Math.min(res[0], act[0]);
				res[1] = Math.max(res[1], act[1]);
			}
		}
		return res;
	}

	public boolean notEmpty(List<Reaction> list) {
		return list.size() != 0;
	}

	public void removeReaction(Reaction re, List<Reaction> removables) {
		inputNetwork.deleteReaction(re);
		removables.remove(re);
	}

	public List<Reaction> removables() {
		List<Reaction> removables = new ArrayList<>();
		List<Reaction> reactions = inputNetwork.getReactionsList();
		for (int i = 0; i < inputNetwork.getReactionsList().size(); i++) {
			Reaction act = reactions.get(i);
			if (!protectedReactions.contains(act)) {
				removables.add(act);
			}
		}
		return removables;
	}

	public Reaction candidateReactForRemoval(List<Reaction> removables) {
		Reaction candidate = null;
		double minVal = Double.MAX_VALUE;
		for (int i = 0; i < removables.size(); i++) {
			Reaction act = removables.get(i);
			double[] fluxRange = calculateFluxRange(act);
			double range = fluxRange[1] - fluxRange[0];
			if (range < minVal) {
				minVal = range;
				candidate = act;
			}
		}
		return candidate;
	}

	public void reinsert(Reaction re) {
		inputNetwork.addReaction(re);
	}

	public int calculateDegreesOfFreedom() {
		return inputNetwork.getReactionsList().size() - rankOfMatrix();
	}

	

	public int rankOfMatrix() {
		// Tomado de: https://www.geeksforgeeks.org/program-for-rank-of-matrix/
		// Author: Utkarsh Trivedi
		int mat[][] = inputNetwork.createAdjacencyMatrix();
		int rank = mat[0].length;
		int R = mat.length;
		for (int row = 0; row < rank; row++) {
			if (mat[row][row] != 0) {
				for (int col = 0; col < R; col++) {
					if (col != row) {
						double mult = (double) mat[col][row] / mat[row][row];

						for (int i = 0; i < rank; i++)

							mat[col][i] -= mult * mat[row][i];
					}
				}
			} else {
				boolean reduce = true;

				for (int i = row + 1; i < R; i++) {
					if (mat[i][row] != 0) {
						swap(mat, row, i, rank);
						reduce = false;
						break;
					}
				}
				if (reduce) {
					rank--;
					for (int i = 0; i < R; i++)
						mat[i][row] = mat[i][rank];
				}
				row--;
			}
		}
		return rank;
	}

	public void swap(int mat[][], int row1, int row2, int col) {
		for (int i = 0; i < col; i++) {
			int temp = mat[row1][i];
			mat[row1][i] = mat[row2][i];
			mat[row2][i] = temp;
		}
	}

	public MetabolicNetwork getInputNetwork() {
		return inputNetwork;
	}

	public void setInputNetwork(MetabolicNetwork inputNetwork) {
		this.inputNetwork = inputNetwork;
	}

	public Set<Metabolite> getProtectedMetabolites() {
		return protectedMetabolites;
	}

	public void setProtectedMetabolites(Set<Metabolite> protectedMetabolites) {
		this.protectedMetabolites = protectedMetabolites;
	}

	public Set<Reaction> getProtectedReactions() {
		return protectedReactions;
	}

	public void setProtectedReactions(Set<Reaction> protectedReactions) {
		this.protectedReactions = protectedReactions;
	}

	public Map<Inequality, Boolean> getProtectedFunctionsAndPhenotypes() {
		return protectedFunctionsAndPhenotypes;
	}

	public void setProtectedFunctionsAndPhenotypes(Map<Inequality, Boolean> protectedFunctionsAndPhenotypes) {
		this.protectedFunctionsAndPhenotypes = protectedFunctionsAndPhenotypes;
	}

	public int getMinimunDegreesOfFreedom() {
		return minimunDegreesOfFreedom;
	}

	public void setMinimunDegreesOfFreedom(int minimunDegreesOfFreedom) {
		this.minimunDegreesOfFreedom = minimunDegreesOfFreedom;
	}

	public int getMinimumNumberOfReactions() {
		return minimumNumberOfReactions;
	}

	public void setMinimumNumberOfReactions(int minimumNumberOfReactions) {
		this.minimumNumberOfReactions = minimumNumberOfReactions;
	}
}
