#include <cudacompare.h>

#include <cstdlib>
#include <vector>


int main()
{
  segment defaultseg;
  defaultseg.raw = std::vector< std::vector<double> >( 50, std::vector<double>(13, 0.5) );
  std::vector<segment> segs( 1300, defaultseg );
  std::vector<std::vector<Match>> matches = findSubSequencesCUDA(segs, true, true);
}
