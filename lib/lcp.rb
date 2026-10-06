# Solves the Linear Complementarity Problem (LCP):
#
#   given M (n x n) and q (n), find z >= 0 and w = Mz + q >= 0
#   such that z_i * w_i == 0 for every i.
#
# Like Simplex, this works directly on an explicit tableau with exact
# Rational arithmetic and exposes every pivot. It runs Lemke's
# complementary pivot algorithm:
#
#   * an artificial variable z0 (with covering vector d = 1) is pivoted
#     into the basis to make the dictionary feasible,
#   * each subsequent entering variable is the complement of the
#     variable that just left, chosen by a lexicographic minimum-ratio
#     test (which also decides degenerate ties and prevents cycling),
#   * the problem is solved exactly when the artificial variable leaves
#     the basis -- never by a floating-point penalty reaching zero,
#   * if the entering column offers no legal pivot, the algorithm stops
#     on a secondary ray and reports the ray as evidence instead of a
#     "solution".
class LCP
  DEFAULT_MAX_PIVOTS = 10_000

  # Index of the artificial variable in the tableau columns.
  ARTIFICIAL = 0

  # The outcome of a solve. +status+ is :solved or :secondary_ray.
  # Iteration records describe every pivot in the complementary path.
  class Result
    attr_reader :status, :z, :w, :iterations, :ray

    def initialize(status:, z: nil, w: nil, iterations: [], ray: nil)
      @status     = status
      @z          = z
      @w          = w
      @iterations = iterations
      @ray        = ray
    end

    def solved?
      @status == :solved
    end
  end

  attr_reader :m, :q, :n, :iterations, :basis
  attr_accessor :max_pivots

  def initialize(m, q)
    @n = q.length
    unless m.length == @n && m.all? {|row| row.length == @n }
      raise ArgumentError, "M must be a square matrix matching the dimension of q"
    end

    @m = m.map {|row| row.map {|v| Rational(v) } }
    @q = q.map {|v| Rational(v) }

    @max_pivots = DEFAULT_MAX_PIVOTS
    setup_tableau
  end

  def self.solve(m, q)
    new(m, q).solve
  end

  # Runs Lemke's algorithm and returns a Result.
  def solve
    setup_tableau
    @iterations = []

    # Trivial case: q >= 0 means z = 0, w = q is already complementary.
    if @rhs.all? {|bi| bi >= 0 }
      return Result.new(:status => :solved,
                        :z => Array.new(@n) { Rational(0) },
                        :w => @q.dup,
                        :iterations => @iterations)
    end

    # Initial pivot: the artificial variable z0 enters the basis in the
    # row with the most negative right-hand side. This single pivot
    # makes the dictionary feasible; no ratio test is involved.
    row     = initial_pivot_row
    leaving = @basis[row]
    do_pivot(row, ARTIFICIAL)
    record_pivot(ARTIFICIAL, leaving, row, nil,
                 "initial pivot: artificial variable z0 enters, " \
                 "#{var_name(leaving)} leaves (most negative right-hand side)")

    # Complementary pivoting: the complement of the variable that just
    # left the basis is the next to enter.
    entering = complement(leaving)
    loop do
      raise "Too many pivots" if @iterations.size >= max_pivots

      row, ratios = ratio_test(entering)
      unless row
        return Result.new(:status => :secondary_ray,
                          :iterations => @iterations,
                          :ray => secondary_ray_evidence(entering))
      end

      leaving = @basis[row]
      do_pivot(row, entering)
      solved  = (leaving == ARTIFICIAL)
      note    = "complementary pivot: #{var_name(entering)} enters " \
                "(complement of #{var_name(complement(entering))}), " \
                "#{var_name(leaving)} leaves by the lexicographic minimum-ratio test"
      note   += "; artificial variable z0 left the basis -- solution found" if solved
      record_pivot(entering, leaving, row, ratios, note)

      if solved
        z, w = extract_solution
        return Result.new(:status => :solved, :z => z, :w => w,
                          :iterations => @iterations)
      end
      entering = complement(leaving)
    end
  end

  # --- tableau -------------------------------------------------------

  # The tableau stores rows of the system
  #
  #   x_B[i] + sum_j T[i][j] * x_j = rhs[i]        (j over all variables)
  #
  # so basic columns are unit columns and the dictionary form is
  # x_B[i] = rhs[i] - sum_j T[i][j] * x_j over the nonbasic j.
  def setup_tableau
    @t = Array.new(@n) { Array.new(num_vars) { Rational(0) } }
    @rhs   = @q.dup
    @basis = (0...@n).map {|i| w_var(i) }
    @n.times do |i|
      @t[i][ARTIFICIAL] = Rational(-1) # covering vector d = 1
      @n.times {|j| @t[i][z_var(j)] = -@m[i][j] }
      @t[i][w_var(i)] = Rational(1)
    end
  end

  def num_vars
    2 * @n + 1
  end

  def z_var(j)
    j + 1
  end

  def w_var(i)
    @n + 1 + i
  end

  def var_name(v)
    if v == ARTIFICIAL
      "z0"
    elsif v <= @n
      "z#{v}"
    else
      "w#{v - @n}"
    end
  end

  # The complement pairs are (z_i, w_i); the artificial variable has none.
  def complement(v)
    return nil if v == ARTIFICIAL
    v <= @n ? w_var(v - 1) : z_var(v - @n - 1)
  end

  def do_pivot(row, entering)
    pivot_value = @t[row][entering]
    new_row     = @t[row].map {|v| v / pivot_value }
    new_rhs     = @rhs[row] / pivot_value

    (0...@n).each do |i|
      next if i == row
      factor = @t[i][entering]
      next if factor == 0
      @rhs[i] -= factor * new_rhs
      num_vars.times {|j| @t[i][j] -= factor * new_row[j] }
    end

    @t[row]     = new_row
    @rhs[row]   = new_rhs
    @basis[row] = entering
  end

  # Row with the most negative rhs; ties are broken lexicographically
  # on the tableau row so the choice is deterministic.
  def initial_pivot_row
    min_rhs    = @rhs.min
    candidates = (0...@n).select {|i| @rhs[i] == min_rhs }
    candidates.min_by {|i| @t[i] }
  end

  # Lexicographic minimum-ratio test over the rows whose entering
  # column entry is positive. Comparing the whole normalized row
  # (rhs first, then every tableau column) keeps the dictionary
  # lexicographically feasible, resolves degenerate ties and
  # guarantees the algorithm cannot cycle.
  def ratio_test(entering)
    candidates = (0...@n).select {|i| @t[i][entering] > 0 }
    return [nil, []] if candidates.empty?

    ratios = candidates.map do |i|
      { :row => i, :basic => var_name(@basis[i]),
        :ratio => @rhs[i] / @t[i][entering] }
    end
    row = candidates.min_by {|i| lexicographic_key(i, entering) }
    [row, ratios]
  end

  def lexicographic_key(row, entering)
    divisor = @t[row][entering]
    ([@rhs[row]] + @t[row]).map {|v| v / divisor }
  end

  def extract_solution
    values = Hash.new(Rational(0))
    @basis.each_with_index {|v, i| values[v] = @rhs[i] }
    z = (0...@n).map {|j| values[z_var(j)] }
    w = (0...@n).map {|i| values[w_var(i)] }
    [z, w]
  end

  # --- reporting -----------------------------------------------------

  def record_pivot(entering, leaving, row, ratios, note)
    @iterations << {
      :step         => @iterations.size + 1,
      :entering     => var_name(entering),
      :leaving      => var_name(leaving),
      :entering_var => entering,
      :leaving_var  => leaving,
      :pivot_row    => row,
      :ratios       => ratios,
      :note         => note,
      :basis        => @basis.map {|v| var_name(v) },
      :basic_values => @rhs.dup,
    }
  end

  def secondary_ray_evidence(entering)
    {
      :entering        => var_name(entering),
      :entering_var    => entering,
      :column          => (0...@n).map {|i| @t[i][entering] },
      :basic_variables => @basis.map {|v| var_name(v) },
      :basic_values    => @rhs.dup,
      :message         => "no positive entry in the column of #{var_name(entering)}: " \
                          "#{var_name(entering)} can grow without bound along a secondary ray, " \
                          "so Lemke's algorithm cannot continue and no solution was found",
    }
  end

  def formatted_tableau
    names  = (0...num_vars).map {|v| var_name(v) }
    values = @t.map {|row| formatted_values(row) }
    rhs    = formatted_values(@rhs)
    width  = (names + values.flatten + rhs + ["basis"]).map(&:size).max
    lines  = [(["basis", "rhs"] + names).map {|s| s.rjust(width) }.join("  ")]
    lines << "-" * lines.first.length
    @n.times do |i|
      lines << ([var_name(@basis[i]), rhs[i]] + values[i]).map {|s| s.rjust(width) }.join("  ")
    end
    lines.join("\n")
  end

  def formatted_values(array)
    array.map {|v| "%2.3f" % v }
  end
end
