/* BridgesSolver.java

   This is a solver for 'Hashiwokakero'
   
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
   exists at that coordinate and -1 if it is an empty space.
   
   An input file can contain an unlimited number of puzzle grids; each will be
   processed separately.
   
   Andrew Stocks - 11/29/2015
*/

import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Scanner;
import java.io.File;

public class BridgesSolver
{
	static final boolean checkEdgeDump = false;
	static final boolean invalidDump = false;
	
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
		Vertex[][] vertices = new Vertex[puzzleGrid.length][puzzleGrid[0].length];
		for (int y = 0; y < puzzleGrid.length; y++)
		{
			for (int x = 0; x < puzzleGrid[0].length; x++)
			{
				if (isVertex(puzzleGrid[y][x]))
				{
					vertices[y][x] = (new Vertex(y,x,puzzleGrid[y][x]));
				}
			}
		}
		
		int solutionEdges = computeCardinality(puzzleGrid)/2; // The number of edges in the solution
		ArrayList<Edge> edgeList = generateEdgeList(puzzleGrid);
		ArrayList<Edge> edges = new ArrayList<Edge>(edgeList);
		ArrayList<Edge> S = new ArrayList<Edge>();
		
		
		// These optimization rules DO NOT properly work when considering nodes that already have edges connected
		// At least for more complicated situations.
		optimizeComputeNeighbours(vertices);
		boolean optimizing = true;
		while (optimizing)
		{
			optimizing = false;
			for (Vertex[] row : vertices)
			{
				for (Vertex v : row)
				{
					if (v == null)
						continue;
					switch (v.currentN)
					{
						// I wrote these in reverse order because my mind was in a different place
						case -1:
						case 0:
							break;
						case 8: // An 8 has all edges used
							optimizeFillEdges(vertices, edges, S, v, v.getNeighbours());
							optimizing = true;
							break;
						case 7: // A 7 as at least one edge in every direction
							optimizeFillEdges(vertices, edges, S, v, v.getUniqueNeighbours());
							optimizing = true;
							break;
						case 6:
							if (v.getDegree() == 6) // Covers cases where a 6 has only 3 neighbours
							{
								optimizeFillEdges(vertices, edges, S, v, v.getNeighbours());
								optimizing = true;
							}
							break;
						case 5:
							if (v.getNeighbourCount() == 3) // A 5 with only 3 neighbours has at least one edge to each
							{
								//System.out.println("Optimizing vertex of degree 5");
								optimizeFillEdges(vertices, edges, S, v, v.getUniqueNeighbours());
								optimizing = true;
							}
							else if (v.getDegree() == 5)
							{
								optimizeFillEdges(vertices, edges, S, v, v.getNeighbours());
								optimizing = true;
							}								
							break;
						case 4:
							if (v.getDegree() == 4) // Exactly two neighbours means two lines to both of them
							{
								//System.out.println("Optimizing vertex of degree 4");
								optimizeFillEdges(vertices, edges, S, v, v.getNeighbours());
								optimizing = true;
							}
							break;
						case 3:
							if (v.getNeighbourCount() == 2) // A 3 with only 2 neighbours will have at least one line to both
							{
								//System.out.println("Optimizing vertex of degree 3");
								optimizeFillEdges(vertices, edges, S, v, v.getUniqueNeighbours());
								optimizing = true;
							}
							else if (v.getDegree() == 3)
							{
								//System.out.println("Optimizing vertex of degree 3");
								optimizeFillEdges(vertices, edges, S, v, v.getNeighbours());
								optimizing = true;
							}
							break;
						case 2:
							if (v.getDegree() == 2)
							{
								//System.out.println("Optimizing vertex of degree 2");
								optimizeFillEdges(vertices, edges, S, v, v.getNeighbours());
								optimizing = true;
							}
							else if (v.n == 2 && v.getNeighbourCount() == 2) // A 2 cannot connect itself to only another 2 or only a 1
							{
								//System.out.println(v + " --- " + v.getUniqueNeighbours());
								ArrayList<Vertex> neighbours = v.getUniqueNeighbours();
								if (neighbours.get(0).n == 1 || neighbours.get(0).n == 2)
								{
									optimizeRemoveEdge(edges, v, neighbours.get(1));
									v.removeNeighbour(neighbours.get(1));
									neighbours.get(1).removeNeighbour(v);
									S.add(optimizeAddEdge(vertices, edges, v, neighbours.get(1)));
								}
								if (neighbours.get(1).n == 1 || neighbours.get(1).n == 2)
								{
									optimizeRemoveEdge(edges, v, neighbours.get(0));
									v.removeNeighbour(neighbours.get(0));
									neighbours.get(0).removeNeighbour(v);
									S.add(optimizeAddEdge(vertices, edges, v, neighbours.get(0)));
								}
							}
							break;
						case 1:
							if (v.getNeighbourCount() == 1)
							{
								//System.out.println("Optimizing vertex of degree 1");
								optimizeFillEdges(vertices, edges, S, v, v.getUniqueNeighbours());
								optimizing = true;
							}
							// else we could check if the other neighbours are all degree 1 except 1
							break;
					}
				}
			}
		}
		
		System.out.println("Optimized: Edges Remaining(" + edges.size() + ") - Removed:" + (edgeList.size() - edges.size()));
		updateGrid(puzzleGrid, edgeList, S);
		dumpGrid(puzzleGrid);
		
		System.out.println("Remaining possible edges:");
		updateGrid(puzzleGrid, edgeList, edges);
		dumpGrid(puzzleGrid);
		RecursiveSolve(puzzleGrid, vertices, edgeList, edges, S, solutionEdges);
		
		return puzzleGrid;
	}
	
	// Find and return the first encountered neighbour
	public static Vertex optimizeFindNeighbour(Vertex[][] vertices, int y, int x, int dY, int dX)
	{
		y += dY;
		x += dX;
		while ((y >= 0 && y < vertices.length && x >= 0 && x < vertices[0].length))
		{			
			//System.out.print(vertices[y][x]);
			if (vertices[y][x] == null)
			{
				y += dY;
				x += dX;
			}
			else if (vertices[y][x].currentN >= 1 && vertices[y][x].currentN <= 8)
			{
				return vertices[y][x];
			}
			else
			{
				break;
			}
		}
		return null;
	}
	
	public static void optimizeFillEdges(Vertex[][] vertices, ArrayList<Edge> edges, ArrayList<Edge> S, Vertex v, ArrayList<Vertex> neighbours)
	{
		for (Vertex neighbour : neighbours)
		{
			optimizeRemoveEdge(edges, v, neighbour);
			v.removeNeighbour(neighbour);
			neighbour.removeNeighbour(v);
			S.add(optimizeAddEdge(vertices, edges, v, neighbour));
		}
	}
	
	// Compute the neighbours for each vertex in the grid
	public static void optimizeComputeNeighbours(Vertex[][] vertices)
	{
		for (int y = 0; y < vertices.length; y++)
		{
			for (int x = 0; x < vertices[0].length; x++)
			{
				
				if (vertices[y][x] != null && vertices[y][x].n >= 1 && vertices[y][x].n <= 8)
				{
					int n = vertices[y][x].n;
					// Check for an edge in each direction
					// Probably more room for submethods in here
					Vertex v = optimizeFindNeighbour(vertices, y, x, 0, -1);
					if (v != null)
					{
						vertices[y][x].addNeighbour(v);
						if (n > 1 && v.n > 1)
							vertices[y][x].addNeighbour(v);
					}
					v = optimizeFindNeighbour(vertices, y, x, 0, 1);
					if (v != null)
					{
						vertices[y][x].addNeighbour(v);
						if (n > 1 && v.n > 1)
							vertices[y][x].addNeighbour(v);
					}
					v = optimizeFindNeighbour(vertices, y, x, -1, 0);
					if (v != null)
					{
						vertices[y][x].addNeighbour(v);
						if (n > 1 && v.n > 1)
							vertices[y][x].addNeighbour(v);
					}
					v = optimizeFindNeighbour(vertices, y, x, 1, 0);
					if (v != null)
					{
						vertices[y][x].addNeighbour(v);
						if (n > 1 && v.n > 1)
							vertices[y][x].addNeighbour(v);
					}
				}
			}
		}
	}
	
	// Will need to remove neighbours that are blocked by the edge
	public static Edge optimizeAddEdge(Vertex[][] vertices, ArrayList<Edge> edges, Vertex start, Vertex end)
	{		
		Edge result = new Edge(start, end, (start.y == end.y) ? 0 : 1);
		
		// For every line of the edge, scan the perpendicular for pairs of vertices
		ArrayList<Vertex> edgeSquares = result.edgeSquares;
		for (int i = 0 ; i < edgeSquares.size(); i++)
		{
			Vertex currentSquare = edgeSquares.get(i);
			// direction == 0 means the edge is vertical, so scan the horizontal
			if (result.direction == 1)
			{
				Vertex left = optimizeFindNeighbour(vertices, currentSquare.y, currentSquare.x, 0, -1);
				Vertex right = optimizeFindNeighbour(vertices, currentSquare.y, currentSquare.x, 0, 1);
				//System.out.println("Left: " + left + " - Right: " + right);
				// If we find a vertex in both directions, they are no longer neighbours of eachother
				if (left != null && right != null)
				{
					left.removeNeighbour(right);
					left.removeNeighbour(right);					
					right.removeNeighbour(left);
					right.removeNeighbour(left);
					optimizeRemoveEdge(edges, left, right);
					optimizeRemoveEdge(edges, left, right);
				}
			}
			else
			{
				Vertex up = optimizeFindNeighbour(vertices, currentSquare.y, currentSquare.x, -1, 0);
				Vertex down = optimizeFindNeighbour(vertices, currentSquare.y, currentSquare.x, 1, 0);
				if (up != null && down != null)
				{
					up.removeNeighbour(down);
					up.removeNeighbour(down);
					down.removeNeighbour(up);
					down.removeNeighbour(up);
					optimizeRemoveEdge(edges, up, down);
					optimizeRemoveEdge(edges, up, down);
				}
			}
		}
		
		start.currentN--;
		end.currentN--;
		
		if (start.currentN == 0)
		{
			for (Vertex neighbour : start.neighbours)
			{
				neighbour.removeNeighbour(start);
				optimizeRemoveEdge(edges, start, neighbour);
			}
		}
		if (end.currentN == 0)
		{
			for (Vertex neighbour : end.neighbours)
			{
				neighbour.removeNeighbour(end);
				optimizeRemoveEdge(edges, end, neighbour);
			}
		}
		
		
		// Decrement the n of both vertices
		// If a vertex reaches 0, remove all its neighbours
		
		// For each square in the edge travel in both directions perpendicular to the edge
		// If a vertex is encountered in both directions (before encountering another edge or the boundary)
		// They are no longer neighbours of eachother.
		return result;
	}
	
	public static void optimizeRemoveEdge(ArrayList<Edge> edgeList, Vertex a, Vertex b)
	{
		for (Edge e : edgeList)
		{
			if ((e.a.y == a.y && e.a.x == a.x && e.b.y == b.y && e.b.x == b.x)
				|| (e.a.y == b.y && e.a.x == b.x && e.b.y == a.y && e.b.x == a.x))
			{
				edgeList.remove(e);
				break;
			}
		}
	}
	
	public static int computeCardinality(int[][] puzzleGrid)
	{
		int cardinality = 0;
		for (int y = 0; y < puzzleGrid.length; y++)
		{
			for (int x = 0; x < puzzleGrid[0].length; x++)
			{
				int square = puzzleGrid[y][x];
				// Found a vertex
				if (isVertex(square))
					cardinality += square;
			}
		}
		return cardinality;
	}
	
	// n is the # of vertices in the puzzle
	public static boolean RecursiveSolve(int[][] puzzleGrid, Vertex[][] vertices, ArrayList<Edge> edgeList, ArrayList<Edge> edges, ArrayList<Edge> S, int n)
	{
		if (S.size() >= n || edges.isEmpty())
		{
			if (S.size() != n)
				return false;
			updateGrid(puzzleGrid, edgeList, S); // Since we're working with an edge list we'll need to make sure the grid we're passing in here is the fully processed one
			/*
			dumpGrid(puzzleGrid);
			for (Edge e: S)
			{
				System.out.println(e);
			}
			*/
			return verifySolution(puzzleGrid);
		}
		ArrayList<Edge> edgesPrime = new ArrayList<Edge>(edges);
		Edge e = edgesPrime.remove(0);
		if (RecursiveSolve(puzzleGrid, vertices, edgeList, edgesPrime, S, n))
			return true;
		// Assume next edge is in the solution:
		S.add(e);
		optimizeAddEdge(vertices, edgesPrime, e.a, e.b);
		if (RecursiveSolve(puzzleGrid, vertices, edgeList, edgesPrime, S, n))
			return true;
		S.remove(e);
		return false;
	}
	
	// Used by generateEdgeList()
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
	
	// Iterate through the puzzleGrid and return a list of all possible edges
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
	public static void updateGrid(int[][] puzzleGrid, ArrayList<Edge> edges, ArrayList<Edge> S)
	{
		// First clear any unused edges
		for (Edge e : edges)
		{
			for (Vertex v : e.edgeSquares)
			{
				if (!isVertex(puzzleGrid[v.y][v.x]))
					puzzleGrid[v.y][v.x] = -1;
			}
		}
		
		for (Edge e : S)
		{
			// Horizontal edge
			if (e.direction == 0)
			{
				for (Vertex v : e.edgeSquares)
				{
					if (puzzleGrid[v.y][v.x] == 10)
						puzzleGrid[v.y][v.x] = 11;
					else if (puzzleGrid[v.y][v.x] == -1)
						puzzleGrid[v.y][v.x] = 10;
				}
			}
			else // Vertical edge
			{
				for (Vertex v : e.edgeSquares)
				{
					if (puzzleGrid[v.y][v.x] == 20)
						puzzleGrid[v.y][v.x] = 21;
					else if (puzzleGrid[v.y][v.x] == -1)
						puzzleGrid[v.y][v.x] = 20;
				}
			}
		}
		
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
				if (!verified[y][x] && puzzleGrid[y][x] > 0)
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
		Used on the puzzleGrid array.
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
	
	public static void dumpGrid(int[][] puzzleGrid)
	{
		for (int[] row : puzzleGrid)
		{
			for (int square : row)
			{
				if (square >= 0 && square <= 8)
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
			System.out.printf("Reading grid %d\n",graphNum);
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
				System.out.printf("Grid for graph %d contains too few values.\n",graphNum);
				System.out.println("Expected " + x*y + " values");
				break;
			}
			System.out.println("Input grid:");
			dumpGrid(G);
			long startTime = System.currentTimeMillis();
			int[][] solvedPuzzle = Solve(G);
			long endTime = System.currentTimeMillis();
			totalTimeSeconds += (endTime-startTime)/1000.0;
						
			if(!verifySolution(solvedPuzzle))
				System.out.println("Invalid solution detected");
			else
			{
				System.out.println("Solution:");
				dumpGrid(solvedPuzzle);
			}
			System.out.println();
		}
		graphNum--;
		System.out.printf("Processed %d grid%s.\n",graphNum,(graphNum != 1)?"s":"");
		System.out.printf("Total Time (seconds): %.2f\n",totalTimeSeconds);
		System.out.printf("Average Time (seconds): %.2f\n",(graphNum>0)?totalTimeSeconds/graphNum:0);
	}
}