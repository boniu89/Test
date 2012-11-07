import java.math.BigDecimal;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;


public class zad {

	
	public static double round(double unrounded, int precision, int roundingMode)
	{
	    BigDecimal bd = new BigDecimal(unrounded);
	    BigDecimal rounded = bd.setScale(precision, roundingMode);
	    return rounded.doubleValue();
	}
	
	//dotproduct
	public static double dp(double[][] a, double[][] b, int z, int x)
	{
		double temp=0;
	
		for(int i=0;i<a.length;i++)
		{
			temp+= a[i][z]*b[i][x];
		}
		return temp;
	}
	
	//dotproduct
	public static double dp(DoubleMatrix2D a, DoubleMatrix2D b, int z, int x)
	{
		double temp=0;
	
		for(int i=0;i<a.rows();i++)
		{
			temp+= a.get(i, z)*b.get(i, x);
		}
		return temp;
	}
	
	//dlugosc
	public static double dl(double[][] a, int z)
	{
		double temp=0;
		
		for(int i=0;i<a.length;i++)
		{
			temp+= a[i][z]*a[i][z];
		}
		double temp2= Math.sqrt(temp);
		return temp2;
	}
	
	//drukowanie macierzy
	public static void printmat(double[][] mat)
	{
		System.out.println("drukowanie macierzy");
		for(int i=0; i < mat.length; i++)
  		{
  			for(int j=0; j<mat[i].length; j++)
  			{
  				System.out.print(mat[i][j] + ", ");
  			}
  			System.out.print("\n");
  		}
	}
	
	//drukowanie macierzy
	public static void printmat(double[] array)
	{
		System.out.println("drukowanie macierzy");
		for(int i=0; i < array.length; i++)
  		{ 			
  			System.out.print(array[i] + ", ");  			
  		}
		System.out.print("\n");
	}
	
	//drukowanie macierzy trojkatnej - koncowej
	public static void printmat2(double[][] mat)
	{
		System.out.println("drukowanie macierzy");
		for(int i=1; i < mat.length; i++)
		{
			for(int j=0; j<mat[i].length-1; j++)
			{
				if(i>j)
				System.out.print(mat[i][j] + ", ");
			}
			System.out.print("\n");
		}
	}
	
	//pod cosinusowe
	public static double sim_cos(double[][] dp, double[] magn , int a, int b)
	{
		double temp=0;
		
		temp = round(dp[a][b]/ (magn[a]*magn[b]),2,BigDecimal.ROUND_HALF_UP);
			
		return temp;
	}
	
	public static double[][] matrix = { 
			{ 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
			{ 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
			{ 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 1.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
			{ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },
			{ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0 },
			{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 },
			{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0 }, 
			};
		
	public static void main(String[] args) {
		// TODO Auto-generated method stub
          DenseDoubleMatrix2D cos = new DenseDoubleMatrix2D(matrix);
          SingularValueDecomposition cos2 = new SingularValueDecomposition(cos);
          DoubleMatrix2D U = cos2.getU();
          DoubleMatrix2D S = cos2.getS();
          DoubleMatrix2D V = cos2.getV();
//          System.out.println("U "+U.toString());
          System.out.println("S "+S.toString());
        // System.out.println("V "+V.toString());
//          
          double[] suma = new double[9];
  		double[] srednia = new double[9];
  		
  		//obliczanie sredniej
  		for(int i=0; i < 9; i++)
  		{
  			suma[i] = 0;
  			srednia[i] = 0;
  			for(int j=0; j<12; j++)
  			{
  				suma[i] += matrix[j][i];
  			}
  			srednia[i] = suma[i] / 12.0;
  			//System.out.print("Œrednia : " + srednia[i] +"\n");
  		}
  		
  	//	wartosc - srednia
  		double[][] matrixNew = new double[12][9];
  		for(int i=0; i < 12; i++)
  		{
  			for(int j=0; j<9; j++)
  			{
  				matrixNew[i][j] = round(matrix[i][j] - srednia[j],2,BigDecimal.ROUND_HALF_UP);
  			//	System.out.print(matrixNew[i][j] + ", ");
  			}
  		//	System.out.print("\n");
  		}
  		
  	  //  obliczanie dotproduct
  	    double[][] dotproduct = new double[9][9];
  	    
  	    for(int i =0;i<9;i++)
  	    {
  	    	for(int j =0;j<9;j++)
  	    	{
  	    	//System.out.println(dp(matrixNew,matrixNew,i,i+1));
  	    		if(i!=j)
  	    		dotproduct[i][j]=dp(matrixNew,matrixNew,i,j);
  	    		//System.out.println(dotproduct[i][j]=dp(matrixNew,matrixNew,i,j));
  	    	}
  	    }
  	    
  	    //obliczanie |A|
  	    double[] magn = new double[9];
  	    
  	    for(int i=0;i<9;i++)
  	    {
  	    	magn[i]=round(Math.sqrt(dp(matrixNew,matrixNew,i,i)),2,BigDecimal.ROUND_HALF_UP);
  	    	
  	    }
  	    
  	    //obliczanie podobienstwa cos
  	    double[][] simcos = new double[9][9];
  	    
  	  for(int i =0;i<9;i++)
	    {
	    	for(int j =0;j<9;j++)
	    	{
	    	//System.out.println(dp(matrixNew,matrixNew,i,i+1));
	    		if(i!=j)
	    		simcos[i][j]=sim_cos(dotproduct, magn,i,j);
	    	}
	    }
  	    System.out.println("Przed SVD");
 	   printmat2(simcos);
 	  
 	   //Obliczanie po SVD
  	//System.out.println("OBLICZANKO2");
  	//System.out.println(V.get(1, 1));
  	  
  	  //Wybranie odpowiednich fragmentow macierzy U S V
  	  DoubleMatrix2D Ux = U.viewPart(0, 0, U.rows(), 2);
  	 // System.out.println(U.toString());
  	//  printmat(Ux);
  	  
  	  DoubleMatrix2D Sx = S.viewPart(0, 0, 2, 2);
  	//  System.out.println(S.toString());
  	 // printmat(Sx);
  	  
  	 DoubleMatrix2D Vx = V.viewPart(0, 0, V.rows(), 2).viewDice();
 	//  System.out.println(V.toString());
 	  
 	//  printmat(Vx);
  	  
  	 
  	 //mnozenie macierzy UxSxV
 	  DoubleMatrix2D Ax = Ux.zMult(Sx, null); 
 	    //System.out.println(Ax);
 	    
 	    Ax= Ax.zMult(Vx,null);
 	    //System.out.println(Ax);
 	   double[][] Axx=Ax.toArray();
 	    
 	   //obliczanie sredniej
 	   double[] suma2 = new double[9];
 		double[] srednia2 = new double[9];
 		
 		for(int i=0; i < 9; i++)
 		{
 			suma2[i] = 0;
 			srednia2[i] = 0;
 			for(int j=0; j<12; j++)
 			{
 				suma2[i] += Axx[j][i];
 			}
 			srednia2[i] = suma2[i] / 12.0;
 			//System.out.print("Œrednia : " + srednia[i] +"\n");
 		}
 		
 	//	obliczanie aartosc - srednia
 		double[][] Axfinal = new double[12][9];
 		for(int i=0; i < 12; i++)
 		{
 			for(int j=0; j<9; j++)
 			{
 				Axfinal[i][j] = /*round(*/Axx[i][j] - srednia2[j];/*,2,BigDecimal.ROUND_HALF_UP);*/
 				//System.out.print(matrixNew[i][j] + ", ");
 			}
 			//System.out.print("\n");
 		}
 	  
 		//obliczanie dotproduct
	    double[][] dotproduct2 = new double[9][9];
	    
	    for(int i =0;i<9;i++)
	    {
	    	for(int j =0;j<9;j++)
	    	{
	    	
	    		if(i!=j)
	    		{
	    		dotproduct2[i][j]=dp(Axfinal,Axfinal,i,j);
	    		//System.out.println(dp(V,V,i,j));
	    		}
	    	}
	    }
	    
	    //obliczanie |A|
	    double[] magn2 = new double[9];
	    
	    for(int i=0;i<9;i++)
	    {
	    	magn2[i]=/*round(*/Math.sqrt(dp(Axfinal,Axfinal,i,i));/*,2,BigDecimal.ROUND_HALF_UP);*/
	    	//System.out.println(magn2[i]);
	    }
	    
	    double[][] simcos2 = new double[9][9];
	    
	  for(int i =0;i<9;i++)
    {
    	for(int j =0;j<9;j++)
    	{
    	//System.out.println(dp(matrixNew,matrixNew,i,i+1));
    		if(i!=j)
    		simcos2[i][j]=sim_cos(dotproduct2, magn2,i,j);
    	}
    }
 	   
	  System.out.println("PO SVD");
 	 //  printmat(dotproduct2);
 	   printmat2(simcos2);
  		  }
			
	}     


