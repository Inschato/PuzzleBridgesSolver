/* BridgesSolver.java

   This is a solver for 'Hashiwokakero'
      
   The input consists of a 2d array puzzle grid in the following format:
       		
   Entry A[i][j] of the grid will be set to 1-8 if a vertex with that number
   exists at that coordinate and -1 if it is an empty space.
   
   It will modify the array in place by adding edges denoted as:
	Single Horizontal		10
	Double Horizontal		11
	Single Vertical			20
	Double Vertical			21
	
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
	
	// Global variables are never wrong, am I right or am I right?
	static int[][] puzzleGrid;
	static Vertex[][] vertices;
	static ArrayList<Edge> edgeList;
	static ArrayList<Edge> S;
	static ArrayList<Vertex> vertexList;
	static int solutionEdges;
	// It's not even using all of these everywhere it could yet.
	
	
	/* Solve(puzzleGrid)
		Given an unsolved puzzle grid, solves it in place
		and returns it.
	 */
	public static int[][] Solve(int[][] G)
	{		
		puzzleGrid = G;
		vertices = new Vertex[puzzleGrid.length][puzzleGrid[0].length];
		vertexList = new ArrayList<Vertex>();
		solutionEdges = 0;
		for (int y = 0; y < puzzleGrid.length; y++)
		{
			for (int x = 0; x < puzzleGrid[0].length; x++)
			{
				if (isValidVertex(puzzleGrid[y][x]))
				{
					vertices[y][x] = (new Vertex(y,x,puzzleGrid[y][x]));
					vertexList.add(vertices[y][x]);
					solutionEdges += puzzleGrid[y][x];
				}
			}
		}		
		solutionEdges /= 2; // We counted each edge twice
		
		// All possible edges for the puzzle
		edgeList = generateEdgeList(vertices);
		// All remaining edges to select from
		ArrayList<Edge> edges = new ArrayList<Edge>(edgeList);
		// Edges in the solution
		S = new ArrayList<Edge>();
		computeNeighbours(vertices);
		
		boolean optimizing = true;
		while (optimizing)
		{
			optimizing = false;	
			for (Vertex v : vertexList)
			{
				if (v.currentN > 0 && v.currentN == v.getAvailableEdgesCount())
				{
					// If a vertex has exactly enough available edges to meet its count, add them
					optimizeFillEdges(vertices, edges, S, v, v.getNeighbours());
					optimizing = true;
				}
				else if (v.currentN > 1 && v.currentN == v.getAvailableEdgesCount()-1)
				{
					// "if the number of slots on me is exactly one less than the sum of my neighbor slots, then add one line to each potential double edge i have" - Hyphenated
					optimizeFillHalfOfDoubleEdges(vertices, edges, S, v);
					optimizing = true;
				}
				else if (v.currentN == 2 && v.n == 2 && v.getNeighbourCount() == 2)
				{
					// A 2 cannot connect itself to only another 2
					ArrayList<Vertex> neighbours = v.getUniqueNeighbours();
					if (neighbours.get(0).n == 2)
					{
						optimizeRemoveEdge(edges, v, neighbours.get(1));
						v.removeNeighbour(neighbours.get(1));
						neighbours.get(1).removeNeighbour(v);
						S.add(optimizeAddEdge(vertices, edges, v, neighbours.get(1)));
						optimizing = true;
					}
					if (neighbours.get(1).n == 2)
					{
						optimizeRemoveEdge(edges, v, neighbours.get(0));
						v.removeNeighbour(neighbours.get(0));
						neighbours.get(0).removeNeighbour(v);
						S.add(optimizeAddEdge(vertices, edges, v, neighbours.get(0)));
						optimizing = true;
					}					
				}
			}
		}
		/*
		if (edges.size() > 1)
		{
			System.out.println("Optimized: Edges Remaining(" + edges.size() + ") - Removed:" + (edgeList.size() - edges.size()));
			updateGrid(S);
			dumpGrid(puzzleGrid);
		
			System.out.println("Remaining possible edges:");
			updateGrid(edges);
			dumpGrid(puzzleGrid);
			
		}
		*/
		// Check if we need to resort to brute forcing the remaining edges
		if (S.size() != solutionEdges)
		{
			RecursiveSolve(edges);		
		}
		updateGrid();
		return puzzleGrid;
	}
	
	// Find and return the first encountered neighbour
	public static Vertex findNeighbour(Vertex[][] vertices, int y, int x, int dY, int dX)
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
	
	// Fill in one of the edges for each side that has two possible edges left
	public static void optimizeFillHalfOfDoubleEdges(Vertex[][] vertices, ArrayList<Edge> edges, ArrayList<Edge> S, Vertex v)
	{
		ArrayList<Vertex> neighbours = v.getNeighbours();		
		// There might be a faster algorithm for finding pairs of the same vertex, this one is n^2 (but n is at most 8)
		for (int i = 0; i < neighbours.size(); i++)
		{
			Vertex a = neighbours.get(i);
			for (int j = i+1; j < neighbours.size(); j++)
			{
				if (a == neighbours.get(j))
				{
					optimizeFillEdges(vertices, edges, S, v, a);
					break;
				}
			}
		}		
	}
	
	public static void optimizeFillEdges(Vertex[][] vertices, ArrayList<Edge> edges, ArrayList<Edge> S, Vertex v, ArrayList<Vertex> neighbours)
	{
		for (Vertex neighbour : neighbours)
		{
			optimizeFillEdges(vertices, edges, S, v, neighbour);			
		}
	}
	
	public static void optimizeFillEdges(Vertex[][] vertices, ArrayList<Edge> edges, ArrayList<Edge> S, Vertex v, Vertex neighbour)
	{
		// Remove the edge from the edge list
		optimizeRemoveEdge(edges, v, neighbour);
		// Remove one instance of the respective vertices from their neighbour lists
		v.removeNeighbour(neighbour);
		neighbour.removeNeighbour(v);
		// Add the edge to the solution list
		S.add(optimizeAddEdge(vertices, edges, v, neighbour));
	}
	
	// Compute the neighbours for each vertex in the grid
	public static void computeNeighbours(Vertex[][] vertices)
	{
		for (int y = 0; y < vertices.length; y++)
		{
			for (int x = 0; x < vertices[0].length; x++)
			{
				
				if (vertices[y][x] != null && vertices[y][x].n >= 1 && vertices[y][x].n <= 8)
				{
					int n = vertices[y][x].n;
					// Check for an edge in each direction
					Vertex[] neighbours = new Vertex[4];
					neighbours[0] = findNeighbour(vertices, y, x, 0, -1);
					neighbours[1] = findNeighbour(vertices, y, x, 0, 1);
					neighbours[2] = findNeighbour(vertices, y, x, -1, 0);
					neighbours[3] = findNeighbour(vertices, y, x, 1, 0);
					for (Vertex neighbour : neighbours)
					{
						if (neighbour != null)
						{
							// If both vertices are 1s, they don't need any edges between them
							if (n > 1 || neighbour.n > 1) 
								vertices[y][x].addNeighbour(neighbour);
							if (n > 1 && neighbour.n > 1)
								vertices[y][x].addNeighbour(neighbour);
						}
					}
				}
			}
		}
	}
	
	public static void optimizeRemoveExtraEdges(ArrayList<Edge> edges, Vertex v)
	{
		ArrayList<Vertex> neighbours = v.getNeighbours();		
		if (v.currentN == 0)
		{
			for (Vertex neighbour : neighbours)
			{
				v.removeNeighbour(neighbour);
				neighbour.removeNeighbour(v);
				optimizeRemoveEdge(edges, v, neighbour);
			}
		}
		else if (v.currentN == 1) // If the vertex only has one edge left to add, turn any remaining possible double edges into singles
		{
			for (int i = 0; i < neighbours.size(); i++)
			{
				Vertex a = neighbours.get(i);
				for (int j = i+1; j < neighbours.size(); j++)
				{
					if (a == neighbours.get(j))
					{
						v.removeNeighbour(a);
						a.removeNeighbour(v);
						optimizeRemoveEdge(edges, v, a);
						break;
					}
				}
			}	
		}
	}
	
	/* optimizeAddEdge()
	 * Adds an edge to the edge list.
	 * Decrements the currentN of each vertex.
	 * Removes edges now blocked by the new edge.
	 * Removes edges no longer possible due to
	 * the reduced n of the vertices.
	 */
	public static Edge optimizeAddEdge(Vertex[][] vertices, ArrayList<Edge> edges, Vertex start, Vertex end)
	{		
		Edge result = new Edge(start, end, (start.y == end.y) ? 0 : 1);
				
		// For each square in the edge travel in both directions perpendicular to the edge
		// If a vertex is encountered in both directions (before encountering another edge or the boundary)
		// They are no longer neighbours of eachother.
		ArrayList<Vertex> edgeSquares = result.edgeSquares;
		for (int i = 0 ; i < edgeSquares.size(); i++)
		{
			Vertex currentSquare = edgeSquares.get(i);
			// direction == 1 means the edge is vertical, so scan the horizontal
			if (result.direction == 1)
			{
				Vertex left = findNeighbour(vertices, currentSquare.y, currentSquare.x, 0, -1);
				Vertex right = findNeighbour(vertices, currentSquare.y, currentSquare.x, 0, 1);
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
				Vertex up = findNeighbour(vertices, currentSquare.y, currentSquare.x, -1, 0);
				Vertex down = findNeighbour(vertices, currentSquare.y, currentSquare.x, 1, 0);
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
		
		// Decrement the n of both vertices
		start.currentN--;
		end.currentN--;
		
		// If a vertex reaches 0, remove all its neighbours
		optimizeRemoveExtraEdges(edges, start);
		optimizeRemoveExtraEdges(edges, end);				
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
	
	private static boolean RecursiveSolve(ArrayList<Edge> edges)
	{
		// Our base case is either when we've run out of edges or when we've got enough in the solution
		// --> This guarentees it will terminate since we're always removing an edge each recursive call
		if (S.size() >= solutionEdges || edges.isEmpty())
		{
			if (S.size() != solutionEdges)
				return false;			
			return verifySolution();
		}
		ArrayList<Edge> edgesPrime = new ArrayList<Edge>(edges);
		Edge e = edgesPrime.remove(0);
		// First assume the next edge is not in the solution
		if (RecursiveSolve(edgesPrime))
			return true;
		// Otherwise assume next edge is in the solution
		S.add(e);
		optimizeAddEdge(vertices, edgesPrime, e.a, e.b);
		if (RecursiveSolve(edgesPrime))
			return true;
		S.remove(e);
		return false;
	}
	
	// Used by generateEdgeList()
	private static Edge findEdge(Vertex[][] vertices, int y, int x, int dY, int dX)
	{
		Vertex currentSquare = vertices[y][x];
		y += dY;
		x += dX;
		while ((y >= 0 && y < vertices.length && x >= 0 && x < vertices[0].length))
		{
			if (vertices[y][x] != null)
			{
				return new Edge(currentSquare, vertices[y][x], (dY != 0) ? 1 : 0);
			}
			y += dY;
			x += dX;
		}
		return null;
	}
	
	// Iterate through the vertices and return a list of all possible edges
	private static ArrayList<Edge> generateEdgeList(Vertex[][] vertices)
	{
		ArrayList<Edge> edgeList = new ArrayList<Edge>();
		for (int y = 0; y < vertices.length; y++)
		{
			for (int x = 0; x < vertices[0].length; x++)
			{
				// Don't add edges outward from a 1 degree vertex
				if (vertices[y][x] != null && vertices[y][x].n > 1)
				{
					// Check for an edge in each direction
					Edge[] edges = new Edge[4];
					edges[0] = findEdge(vertices, y, x, 0, -1);
					edges[1] = findEdge(vertices, y, x, 0, 1);
					edges[2] = findEdge(vertices, y, x, -1, 0);
					edges[3] = findEdge(vertices, y, x, 1, 0);
					for (Edge e : edges)
					{
						if (e != null)
							edgeList.add(e);
					}
				}
			}
		}
		return edgeList;
	}			

	/* verifySolution()
	 * Checks that the current solution S is valid by doing a
	 * traversal to ensure all vertices are connected and checking
	 * that each one has the correct number of adjacent edges
	 */
	private static boolean verifySolution()
	{		
		ArrayList<Vertex> queue = new ArrayList<Vertex>();
		ArrayList<Vertex> visited = new ArrayList<Vertex>();
		queue.add(vertexList.get(0));
		visited.add(vertexList.get(0));
		while (!queue.isEmpty())
		{			
			Vertex currentVertex = queue.remove(0);
			int n = 0;			
			
			for (Edge e : S)
			{				
				if (currentVertex == e.a)
				{
					n++;
					if (!visited.contains(e.b))
					{
						queue.add(e.b);
						visited.add(e.b);						
					}
				}
				else if (currentVertex == e.b)
				{
					n++;
					if (!visited.contains(e.a))
					{
						queue.add(e.a);
						visited.add(e.a);					
					}
				}
			}						
			if (currentVertex.n != n)
				return false;
		}		
		return visited.size() == vertexList.size();
	}
	
	// Returns true if n is between 1 and 8, inclusive.
	// These are the possible values for a vertex (island)
	private static boolean isValidVertex(int n)
	{
		return (n >= 1 && n <= 8);
	}
	
	// Update the grid by removing edges in 'edges' and adding edges in S	
	// This is the final step before returning the solution.
	private static void updateGrid()
	{
		// First clear any unused edges
		for (Edge e : edgeList)
		{
			for (Vertex v : e.edgeSquares)
			{
				if (!isValidVertex(puzzleGrid[v.y][v.x]))
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
	
	// Prints an ACSII representation of the puzzle to the console
	private static void dumpGrid(int[][] puzzleGrid)
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
}