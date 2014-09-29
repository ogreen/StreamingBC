StreamingBC
===========


This repository contains the algorithm I designed for computing betweenness centrality for dynamic graphs.
For this short README file, I will not go into the algorithmic details; rather I will focus on details relevant for any developer using this code.

FIRST OF ALL, under no warranty is this code BUG-FREE. This code was used mostly by me for developing the original dynamic graph algorithm and then extending it to support both parallelism and approximation. As I used it for testing purposes, I would comment in/out whole sections of code to test different hypothesises. In all likelihood this is when bugs were introduced.

The original storage bounds presented in Brandes's 2001 seminal paper has an upperbound of O(V+E) memory. The work on the dynamic graph algorithm led me to believe that this is bound is to high and can be improved by making a simple modifcation to Brandes's algorithm - continue using the parent-child approach without storing the parents. Meaning, a child must search all of its neighbors and figure out which of these are its parent - this is a simply query on the distance to the root for each of the adjacent neighbors. This insight led to a reduced storage complexity of O(V) and allows for better parallelism.

Currently the main file supports graphs from the DIMACS 10 challenge. It reads these graphs from file to memory and stores them in CSR format. Following this, the CSR graph is copied into a STINGER graph instance. STINGER is a graph data structure developed at GT by some awesome people including Rob McColl and David Ediger. STINGER allows for fast graph updates without the constrains of CSR. This repository uses an old version of STINGER.
For a newer version of STINGER go to: https://github.com/robmccoll/stinger
 

Once the graph is in STINGER, that is where the "magic" happens. For those of you interested in using CSR, updating the code to support CSR instead of STINGER requires little work - simply change the STINGER loop directives.

The executable supports multiple command-line parameters:
-N file_name
-T thread_count
-K number_roots_for_approximation  (When K==|V|, the exact BC is computed).
-R seed  ( This controls the seed for the random number generator - useful when you need deterministic results. When seed==0, the time is used as the seed.)

-O operation_type (At different points in time, the code supported both insertions and deletions. Currently, the code only supports insertion. In recent times I have discovered that the deletions and insertions are very similar. Use operation_type=0.


In the code you will see support for additional parameters, in all likelihood, these are no longer in use or were part of other main tester functions that I had at one point or another.

 
I have added below the usual disclaimer just to make it clear that I am not responsible in any way for the behavior of any part of the code (though I understand/support/proved that the algorithm is correct):

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
