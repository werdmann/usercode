#include <limits>
#include <queue>


class FFA{
  // Ford-Fulkerson-Algorithm / Edmond-Karp version, courtesy of chatgpt
public:
  int n;
  unsigned int nrec, nsim;
  std::vector<std::vector<int>> capacity, flow;
  std::vector<int> parent, dist;

  int breadth_first_search(int s, int t) {
    fill(parent.begin(), parent.end(), -1);
    fill(dist.begin(), dist.end(), std::numeric_limits<int>::max());
    std::queue<int> q;
    q.push(s);
    parent[s] = -2;
    dist[s] = 0;
    while (!q.empty()) {
      int u = q.front();
      q.pop();
      for (int v = 0; v < n; v++) {
	if (parent[v] == -1 && capacity[u][v] - flow[u][v] > 0) {
	  parent[v] = u;
	  dist[v] = dist[u] + 1;
	  q.push(v);
	}
      }
    }
    return parent[t] != -1;
  }
  
  int depth_first_search(int u, int t, int f) {
    if (u == t) {
      return f;
    }
    for (int v = 0; v < n; v++) {
      if (dist[v] == dist[u] + 1 && capacity[u][v] - flow[u][v] > 0) {
	int df = depth_first_search(v, t, std::min(f, capacity[u][v] - flow[u][v]));
	if (df > 0) {
	  flow[u][v] += df;
	  flow[v][u] -= df;
	  return df;
	}
      }
    }
    return 0;
  }
  
  int maxflow(int s, int t){
    int totalFlow = 0;
    while (breadth_first_search(s, t)) {
      int df;
      while ((df = depth_first_search(s, t,  std::numeric_limits<int>::max()))) {
	totalFlow += df;
      }
    }
    return totalFlow;
  }


  FFA(unsigned int nrec, unsigned int nsim, std::vector<std::pair<unsigned int,unsigned int>> links, int & nmatch): nrec(nrec), nsim(nsim)
  {
    n = nrec + nsim + 2; // add source (0) and sink  (n-1)
    capacity.assign(n, std::vector<int>(n, 0));
    for(unsigned int i =  0; i < nrec; i++) capacity[0][1+i] = 1;
    for(unsigned int i =  0; i < nsim; i++) capacity[1+nrec+i][n-1] = 1;
    for(auto recsim :  links){
      unsigned int rec = recsim.first;
      unsigned int sim = recsim.second;
      capacity[1+rec][1+nrec+sim] = 1;
    }

    // Set up the flow network.
    flow.assign(n, std::vector<int>(n, 0));
    parent.resize(n);
    dist.resize(n);

    nmatch = maxflow(0, n-1);
  }

  std::vector<std::pair<unsigned int,unsigned int>> get_assignment(){ 
    // return an assignment, i.e. matched pairs (rec, sim), 
    // note: the assignment returned is not necessarily unique and often not the "best" 
    //       it is just one possible matching that results in the maximum flow
    //       this algorithm has no concept of "cost", all allowed links are considered equally good
    std::vector<std::pair<unsigned int,unsigned int>> matches;
    for (unsigned int k = 0; k < nrec; k++){
      for (unsigned int j = 0; j < nsim; j++) {
	auto f = flow[k + 1][j + 1 + nrec];
	if (f > 0.99){
	  matches.emplace_back(k,j);
	}
      }
    }
    return matches;
  }

};


/*
int main() {

    int nrec = 4, nsim = 5;
    std::vector<std::pair<unsigned int,unsigned int>> a;
    a.emplace_back(0,0);
    a.emplace_back(0,1);
    a.emplace_back(1,1);
    a.emplace_back(0,2);
    a.emplace_back(2,3);
    a.emplace_back(3,3);

    int nmatch;
    FFA match = FFA(nrec, nsim, a, nmatch);
    std::cout << "Maximum flow: " << nmatch << std::endl;
    auto assignment = match.get_assignment();
    std::cout << "test  " << (nmatch == assignment.size()) << std::endl;
    for(auto & recsim : assignment){
      std::cout << recsim.first << " : " << recsim.second << std::endl;
    }

    return 0;
}
*/
