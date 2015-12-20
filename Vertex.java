import java.util.ArrayList;
import java.util.HashSet;

public class Vertex
{
	public ArrayList<Vertex> neighbours;
	public HashSet<Vertex> uniqueNeighbours;
	public int y;
	public int x;
	public int n; // The original number of this vertex
	public int currentN; // The modified number of this vertex
	
	public Vertex(int y, int x, int n)
	{
		neighbours = new ArrayList<Vertex>();
		uniqueNeighbours = new HashSet<Vertex>();
		this.y = y;
		this.x = x;			
		this.n = n;
		currentN = n;
	}
	
	public void addNeighbour(Vertex v)
	{
		neighbours.add(v);
		uniqueNeighbours.add(v);
	}
	
	public void removeNeighbour(Vertex v)
	{
		neighbours.remove(v);
		if (!neighbours.contains(v))
		{
			uniqueNeighbours.remove(v);
		}
	}					
	
	public ArrayList<Vertex> getNeighbours()
	{
		return new ArrayList<Vertex>(neighbours);
	}
	
	public ArrayList<Vertex> getUniqueNeighbours()
	{
		return new ArrayList<Vertex>(uniqueNeighbours);
	}
	
	public int getDegree()
	{
		return neighbours.size();
	}
		
	public int getNeighbourCount()
	{
		return uniqueNeighbours.size();
	}
	
	public String toString()
	{
		return "V: " + x + "," + y + " - " + n;
	}
}