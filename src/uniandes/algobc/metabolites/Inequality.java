package uniandes.algobc.metabolites;

public class Inequality {
	private Double LHSConstant;
	private Double RHSConstant;
	private Double[] LHSValues;
	private Double[] RHSValues;

	public Inequality() {
		LHSConstant = new Double(1);
		RHSConstant = new Double(1);
	}

	public void createInequality(Double LHSConstant, Double RHSConstant, Double delta, Double lowerLimitLHS,
			Double upperLimitLHS, Double lowerLimitRHS, Double upperLimitRHS) {
		this.LHSConstant = LHSConstant;
		this.RHSConstant = RHSConstant;
		this.LHSValues = new Double[(int) Math.ceil((upperLimitLHS - lowerLimitLHS) / delta)];
		this.RHSValues = new Double[(int) Math.ceil((upperLimitRHS - lowerLimitRHS) / delta)];
		double lhsAcum=lowerLimitLHS;
		double rhsAcum = lowerLimitRHS;
		for (int i = 0; i < LHSValues.length; i++) {
			LHSValues[i] = lhsAcum + delta;
		}
		for (int i = 0; i < RHSValues.length; i++) {
			RHSValues[i] = rhsAcum + delta;
		}
	}
	
	public void createInequality(Double LHSConstant, Double RHSConstant, Double delta, Double lowerLimitLHS,Double upperLimitLHS ) {
		createInequality(LHSConstant, RHSConstant, delta, lowerLimitLHS,upperLimitLHS,  new Double(0),new Double(0));
	}
	
	public boolean checkFeasibility() {
		boolean b = true;
		if(RHSValues.length >0) {
			for (int i = 0; i < LHSValues.length; i++) {
				b = b && LHSConstant*LHSValues[i] <=RHSConstant*RHSValues[i];
			}
		}
		else {
			for (int i = 0; i < LHSValues.length; i++) {
				b = b && LHSConstant*LHSValues[i] <=RHSConstant;
			}
		}
		return b; 

	}
	public double[] getFeasibilityRange(double lowerLimit, double upperLimit, double delta) {
		double[] res = new double[2];
		double[] values = new double[(int) Math.ceil((upperLimit - lowerLimit) / delta)];
		for (int i = 0; i < LHSValues.length; i++) {
			values[i] = lowerLimit + delta;
		}
		int i =0;
		if(RHSValues.length >0) {
			while(!(LHSConstant*values[i] <=RHSConstant*values[i])) {
				i++;
			}
			res[0]=values[i];
			i = values.length-1;
			while(values[i]>res[0] && !(LHSConstant*values[i] <=RHSConstant*values[i])) {
				i--;
			}
			res[1] = values[i];
		}
		else {
			while(!(LHSConstant*values[i] <=RHSConstant)) {
				i++;
			}
			res[0]=values[i];
			i = values.length-1;
			while(values[i]>res[0]&&!(LHSConstant*values[i] <=RHSConstant)) {
				i--;
			}
			res[1] = values[i];
		}
		return res;
	}
	public double getLHSConstant() {
		return LHSConstant;
	}

	public void setLHSConstant(double lHSConstant) {
		LHSConstant = lHSConstant;
	}

	public double getRHSConstant() {
		return RHSConstant;
	}

	public void setRHSConstant(double rHSConstant) {
		RHSConstant = rHSConstant;
	}

	public Double[] getLHSValues() {
		return LHSValues;
	}

	public void setLHSValues(Double[] lHSValues) {
		LHSValues = lHSValues;
	}

	public Double[] getRHSValues() {
		return RHSValues;
	}

	public void setRHSValues(Double[] rHSValues) {
		RHSValues = rHSValues;
	}
}
