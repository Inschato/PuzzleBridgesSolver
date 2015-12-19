/* BridgesSolver.java   
   To interactively provide test inputs, run the program with
	java BridgesSolver
	
   To conveniently test the algorithm with a large input, create a text file
   containing one or more test puzzle grids (in the format described below) and run
   the program with
	java BridgesSolver file.txt
   where file.txt is replaced by the name of the text file.
   
   The input consists of a series of puzzle grids in the following format:
   
    <number of columns> <number of rows>
	<row 1>
	...
	<row n>
	
   Entry A[i][j] of the grid will be set to 1-8 if a vertex with that number
   exists at that coordinate and . if it is an empty space.
   
   An input file can contain an unlimited number of puzzle grids; each will be
   processed separately.
   
   Andrew Stocks - 11/29/2015
*/

import java.util.Arrays;
import java.util.HashSet;
import java.util.ArrayList;
import java.util.Scanner;
import java.io.File;

public class BridgesSolver
{
	static boolean checkEdgeDump = false;
	static boolean invalidDump = false;
	
	/* Solve(puzzleGrid)
		Given an unsolved puzzle grid, returns a solved puzzle grid.
		Edges are denoted with the following values:
		Single Horizontal		10
		Double Horizontal		11
		Single Vertical			20
		Double Vertical			21
	 */
	public static int[][] Solve(int[][] puzzleGrid)
	{
		dumpGrid(puzzleGrid);
		ArrayList<Edge> edgeList = generateEdgeList(puzzleGrid);
		for (Edge e: edgeList)
		{
			System.out.println(e);
		}
		return puzzleGrid;
	}
	
	// n is the # of vertices in the puzzle
	public static boolean RecursiveSolve(int[][] puzzleGrid, ArrayList<Edge> edges, HashSet<Edge> S, int n)
	{
		if (S.size() >= n) 
		{
			updateGrid(puzzleGrid, edges, S); // Since we're working with an edge list we'll need to make sure the grid we're passing in here is the fully processed one
			return verifySolution(puzzleGrid);
		}
		return false;
	}
	
	public static Edge findEdge(int[][] puzzleGrid, int y, int x, int dY, int dX)
	{
		Vertex currentSquare = new Vertex(y, x, puzzleGrid[y][x]);
		y += dY;
		x += dX;		
		while (inBounds(puzzleGrid, y, x))
		{
			if (isVertex(puzzleGrid[y][x]))
			{
				return new Edge(currentSquare, new Vertex(y, x, puzzleGrid[y][x]), (dY != 0) ? 1 : 0);
			}
			y += dY;
			x += dX;
		}
		return null;
	}
	
	public static ArrayList<Edge> generateEdgeList(int[][] puzzleGrid)
	{
		ArrayList<Edge> edgeList = new ArrayList<Edge>();
		for (int y = 0; y < puzzleGrid.length; y++)
		{
			for (int x = 0; x < puzzleGrid[0].length; x++)
			{
				if (puzzleGrid[y][x] >= 2 && puzzleGrid[y][x] <= 8)
				{					
					// Check for an edge in each direction					
					Edge e = findEdge(puzzleGrid, y, x, 0, -1);
					if (e != null)
						edgeList.add(e);
					e = findEdge(puzzleGrid, y, x, 0, 1);
					if (e != null)
						edgeList.add(e);
					e = findEdge(puzzleGrid, y, x, -1, 0);
					if (e != null)
						edgeList.add(e);
					e = findEdge(puzzleGrid, y, x, 1, 0);
					if (e != null)
						edgeList.add(e);
				}
			}
		}
		return edgeList;
	}
	
	// Update the grid by removing edges in 'edges' and adding edges in S
	// TODO: Figure out how to differentiate between 1-2 edges in the same spot (Or is this just handled when we create the edge list)
	public static void updateGrid(int[][] puzzleGrid, ArrayList<Edge> edges, HashSet<Edge> S)
	{
		
	}

	public static boolean verifySolution(int[][] puzzleGrid)
	{
		// Use this so after we've checked every connected vertex we check that every other space is empty.
		boolean[][] verified = new boolean[puzzleGrid.length][puzzleGrid[0].length];
		boolean foundStart = false;
		ArrayList<Vertex> queue = new ArrayList<Vertex>();
		
		// Just find a starting vertex and then do a breadth first search
		// This way we also verify that everything is connected
		for (int y = 0; y < puzzleGrid.length; y++)
		{
			for (int x = 0; x < puzzleGrid[0].length; x++)
			{				
				int square = puzzleGrid[y][x];
				// Found a vertex
				if (isVertex(square))
				{
					verified[y][x] = true;
					queue.add(new Vertex(y, x, square));
					foundStart = true;
					break;
				}
			}
			if (foundStart)
				break;
		}

		Vertex currentVertex;
		while (!queue.isEmpty())
		{
			currentVertex = queue.remove(0);
			int y = currentVertex.y;
			int x = currentVertex.x;
			//System.out.println("Checking Vertex at " + y + ", " + x);
			int leftWeight = verifyAdjacent(puzzleGrid, verified, queue, y, x, 0, -1);
			int rightWeight = verifyAdjacent(puzzleGrid, verified, queue, y, x, 0, 1);
			int upWeight = verifyAdjacent(puzzleGrid, verified, queue, y, x, -1, 0);
			int downWeight = verifyAdjacent(puzzleGrid, verified, queue, y, x, 1, 0);
			
			if (currentVertex.n == leftWeight + rightWeight + upWeight + downWeight)
				verified[y][x] = true;
			else
			{			
				if (invalidDump)
				{
					System.out.println("Vertex at " + x + ", " + y + " has invalid edges");
					System.out.println("Degree: " + currentVertex.n + " , Edges: " + leftWeight + ", " + rightWeight + ", " + upWeight + ", " + downWeight);
				}
				return false;
			}
		}
		// Every unverified square should be an empty space.
		for (int y = 0; y < puzzleGrid.length; y++)
		{
			for (int x = 0; x < puzzleGrid[0].length; x++)
			{
				if (!verified[y][x] && puzzleGrid[y][x] != -1)
				{
					if(invalidDump)
						System.out.println("Invalid square at " + x + ", " + y + " --> " + puzzleGrid[y][x]);
					return false;
				}
			}
		}
		return true;
	}
	
	/* verifyEdge()
		Given an edge type and direction continue along the direction until reaching a vertex.
		If it encounters anything else (such as a different kind of edge
		or the boundary of the grid) it returns false.
	 */
	public static boolean verifyEdge(int[][] puzzleGrid, boolean[][] verified, ArrayList<Vertex> queue, int y, int x, int dY, int dX, int edgeType)
	{		
		while (inBounds(puzzleGrid, y, x))
		{
			int currentSquare = puzzleGrid[y][x];
			if (isVertex(currentSquare))
			{
				// If we encounter a vertex for the first time add it to the BFS queue
				if (!verified[y][x])
				{
					verified[y][x] = true;
					queue.add(new Vertex(y, x, puzzleGrid[y][x]));
				}
				return true;
			}
			else if (currentSquare != edgeType)
			{
				break;
			}
			verified[y][x] = true;
			y += dY;
			x += dX;			
		}
		return false;
	}
		
	// Assumption: No vertices in directly adjacent squares.
	public static int verifyAdjacent(int[][] puzzleGrid, boolean[][] verified, ArrayList<Vertex> queue, int y, int x, int dY, int dX)
	{
		if (checkEdgeDump)
			System.out.println("Checking edge");
		int weight = 0;
		y += dY;
		x += dX;
		if (inBounds(puzzleGrid, y, x))
		{
			int currentSquare = puzzleGrid[y][x];
			if (((currentSquare == 10 || currentSquare == 11) && (dX == 1 || dX == -1))
				|| ((currentSquare == 20 || currentSquare == 21) && (dY == 1 || dY == -1)))
			{
				if (verifyEdge(puzzleGrid, verified, queue, y, x, dY, dX, currentSquare))
					weight = currentSquare % 10 + 1;
			}
		}
		return weight;
	}
		
	// Returns true if the y and x indexes are inside the array
	public static boolean inBounds(int[][] puzzleGrid, int y, int x)
	{
		return (y >= 0 && y < puzzleGrid.length && x >= 0 && x < puzzleGrid[0].length);
	}	
		
	public static boolean isVertex(int n)
	{
		return (n >= 1 && n <= 8);
	}		
	
	public static class Vertex
	{
		public ArrayList<Vertex> neighbours;
		public int y;
		public int x;		
		public int n;
		
		public Vertex(int y, int x, int n)
		{
			neighbours = new ArrayList<Vertex>();
			this.y = y;
			this.x = x;			
			this.n = n;
		}
		
		public void addNeighbour(Vertex v)
		{
			neighbours.add(v);			
		}
		
		public void removeNeighbour(Vertex v)
		{
			neighbours.remove(v);
		}					
		
		public int getDegree()
		{
			return neighbours.size();
		}
	}
	
	public static class Edge
	{
		public Vertex a;
		public Vertex b;
		public int direction;
		public ArrayList<Vertex> edgeSquares;
		
		public Edge(Vertex a, Vertex b, int direction)
		{
			this.a = a;
			this.b = b;
			this.direction = direction;
			edgeSquares = new ArrayList<Vertex>();
			addLine(a.y, a.x, b.y, b.x);
			removeSquare(a);
			removeSquare(b);
		}
		
		public void addSquare(Vertex v)
		{
			edgeSquares.add(v);
		}
		
		public void removeSquare(Vertex v)
		{
			edgeSquares.remove(v);
		}
		
		public void addLine(int y1, int x1, int y2, int x2)
		{
			if (y1 == y2)
			{
				while (x1 > x2)
				{
					addSquare(new Vertex(y1, x1--, 0));
				}
				while (x1 < x2)
				{
					addSquare(new Vertex(y1, x1++, 0));
				}
			}
			else
			{
				while (y1 > y2)
				{
					addSquare(new Vertex(y1--, x1, 0));
				}
				while (y1 < y2)
				{
					addSquare(new Vertex(y1++, x1, 0));
				}
			}			
		}
		
		public String toString()
		{
			return "E: " + a.x + "," + a.y + " - " + b.x + "," + b.y;
		}
	}
	
	public static void dumpGrid(int[][] puzzleGrid)
	{
		for (int[] row : puzzleGrid)
		{
			for (int square : row)
			{
				if (square >= 1 && square <= 8)
					System.out.print(square);
				else if (square == -1)
					System.out.print(' ');
				else if (square == 10)
					System.out.print("-");
				else if (square == 11)
					System.out.print("=");
				else if (square == 20)
					System.out.print("|");
				else if (square == 21)
					System.out.print("H");
				else
					System.out.print('?');
			}
			System.out.println();
		}
	}
	
	// Input code adapted from code by Bill Bird (~2015, UVic)
	public static void main(String[] args)
	{
		Scanner s;
		if (args.length > 0)
		{
			try
			{
				s = new Scanner(new File(args[0]));
			}
			catch(java.io.FileNotFoundException e)
			{
				System.out.printf("Unable to open %s\n",args[0]);
				return;
			}
			System.out.printf("Reading input values from %s.\n",args[0]);
		}
		else
		{
			s = new Scanner(System.in);
			System.out.printf("Reading input values from stdin.\n");
		}
		
		int graphNum = 0;
		double totalTimeSeconds = 0;
		
		//Read graphs until EOF is encountered (or an error occurs)
		while(true)
		{
			graphNum++;
			if(graphNum != 1 && !s.hasNextInt())
				break;
			System.out.printf("Reading graph %d\n",graphNum);
			int x = s.nextInt();
			int y = s.nextInt();
			int[][] G = new int[y][x];
			int valuesRead = 0;
			for (int i = 0; i < y && s.hasNextInt(); i++)
			{
				for (int j = 0; j < x && s.hasNextInt(); j++)
				{
					G[i][j] = s.nextInt();
					valuesRead++;
				}
			}
			if (valuesRead < x*y)
			{
				System.out.printf("Adjacency matrix for graph %d contains too few values.\n",graphNum);
				System.out.println("Expected " + x*y + " values");				
				break;
			}
			long startTime = System.currentTimeMillis();			
			int[][] solvedPuzzle = Solve(G);
			long endTime = System.currentTimeMillis();
			totalTimeSeconds += (endTime-startTime)/1000.0;
						
			if(!verifySolution(solvedPuzzle))
				System.out.println("Invalid Solution Detected");
			System.out.println();
		}
		graphNum--;
		System.out.printf("Processed %d graph%s.\n",graphNum,(graphNum != 1)?"s":"");
		System.out.printf("Total Time (seconds): %.2f\n",totalTimeSeconds);
		System.out.printf("Average Time (seconds): %.2f\n",(graphNum>0)?totalTimeSeconds/graphNum:0);		
	}
}