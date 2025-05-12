using ModelingToolkit

Start_time = 0.0
Time_step = 0.01
Stop_time = 30.0
tspan = (Start_time, Stop_time)

@parameters t
D = Differential(t)

@connector Pin begin
  p(t)
  q(t), [connect = Flow]
  T(t), [connect = Stream]
end

@mtkmodel OnePort begin
  @components begin
    out = Pin()
    in = Pin()
  end
  @variables begin
    Δp(t)
    q(t)
    T(t)
  end
  @equations begin
    Δp ~ out.p - in.p
    0 ~ in.q + out.q
    q ~ in.q
    in.T ~ instream(in.T)
    out.T ~ T
    D(T) ~ q * (T - in.T)
  end
end

@mtkmodel Resistor begin
  @extend OnePort()

  @parameters begin
    R = 1.0
  end
  @equations begin
      Δp ~ -q * R
  end
end

@mtkmodel Source begin
  @parameters begin
    T_source = 100.0
  end
  @components begin
    out = Pin()
  end
  @equations begin
    # out.q ~ 1.0
    out.p ~ 10.0
    out.T ~ T_source
  end
end


@mtkmodel Junction begin
  @components begin
    in = Pin()
    out1 = Pin()
    out2 = Pin()
  end
  @equations begin
    in.p ~ out1.p
    in.p ~ out2.p
    0 ~ in.q + out1.q + out2.q
    out1.T ~ instream(in.T)
    out2.T ~ instream(in.T)
  end
end

@named source = Source()
@named junction = Junction()
@named r1 = Resistor(R=1.0)
@named r2 = Resistor(R=2.0)

eqs = [
  # connect(source.out, junction.in),
  # connect(junction.out1, r1.in),
  # connect(junction.out2, r2.in)
  connect(source.out, r1.in, r2.in)
]

@named model = ODESystem(eqs, t)
connected = compose(model, [source, junction, r1, r2])
expanded = expand_connections(connected)

for u in equations(expand_connections(connected)) # Debuggine
    println(u)
end

circ_sys = structural_simplify(connected)

#### Debugging
# These 3 lines check the total number of equations and unknowns in the simplified system, as well as the number of equations in the original system.
equations(expand(circ_sys))
unknowns(circ_sys)

prob = ODEProblem(circ_sys, [], tspan)

@time Sol = solve(prob, Tsit5(), reltol=1e-10, abstol=1e-12, maxiters=1e8, saveat=Start_time:Time_step:Stop_time)

display(Plot)