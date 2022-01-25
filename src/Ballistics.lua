--[[
    class BallisticUtilities

    The class that contains the internal processing functions to produce the right numbers for the BallisticsFunctions class. This can be accessed through BallisticsFunctions via BallisticsFunctions.Utilities.
]]
--[=[
    @prop Precision number
    @within BallisticsFunctions
    @private
    The target precision used to compare numbers and measure simularity in [BallisticsFunctions:CompareNumbers]. Set to 1e-3 (or 0.001) by default.
]=]
--[=[
    @prop MaxNumber number
    @within BallisticsFunctions
    @private
    The MaxNumber used to replace infinity (or math.huge) during [BallisticsFunctions:ProduceEstimate]. Set to 2147483646 by default.
]=]
local BallisticUtilities = {
    Precision = 1e-3;
    MaxNumber = 2147483646;
} do
    BallisticUtilities.__index = BallisticUtilities
    BallisticUtilities = setmetatable(BallisticUtilities,BallisticUtilities)
end
--[=[
    @within BallisticsFunctions
    @private
    A function that returns the values of coefficients required to identify the polynomials of an equation
    and find the solutions to said polynomials.
    
    Process found here: https://docs.google.com/document/d/1TKhiXzLMHVjDPX3a3U0uMvaiW1jWQWUmYpICjIDeMSA/edit
    
    @param ProjectileSpeed -- The speed that the projectile will have when fired.
    @param DeltaPosition -- The overall position where the shooter's at (0,0,0).
    @param DeltaVelocity -- The overall velocity where the shooter's at (0,0,0).
    @param DeltaAcceleration -- The overall acceleration where the shooter's at (0,0,0).
    @param DeltaJerk -- The overall jerk where the shooter's at (0,0,0).
    @return number, number?, number?, number?, number?, number?, number? -- Returns the coefficients in the format of T0, T1, T2, T3, T4, T5, T6, where: T0 + T1x + T2x^2 + T3x^3 + T4x^4 + T5x^5 + T6x^6 = 0
]=]
function BallisticUtilities:ProduceCoefficients(ProjectileSpeed: number, DeltaPosition:Vector3, DeltaVelocity:Vector3?, DeltaAcceleration:Vector3?, DeltaJerk:Vector3?): (number, number|nil, number|nil, number|nil, number|nil, number|nil, number|nil) --https://docs.google.com/document/d/1TKhiXzLMHVjDPX3a3U0uMvaiW1jWQWUmYpICjIDeMSA/edit
    local T0: number do
        local PSquare = DeltaPosition*DeltaPosition
        T0 = (PSquare.X+PSquare.Y+PSquare.Z)
    end
    if DeltaVelocity ~= nil then
        local T1: number do
            T1 = 2*DeltaPosition:Dot(DeltaVelocity)
        end
        local T2: number do --(Ax^2 + Ay^2 + Az^2)/4
            local VSquare = DeltaVelocity*DeltaVelocity
            T2 = VSquare.X + VSquare.Y + VSquare.Z - ProjectileSpeed*ProjectileSpeed
        end
        if DeltaAcceleration == nil then
            return T0,T1,T2
        else
            do --T2 = (Vx^2 + Px*Ax + Vy^2 + Py*Ay + Vz^2 + Pz*Az - s^2)
                local PTimesA = DeltaPosition*DeltaAcceleration
                T2 += PTimesA.X + PTimesA.Y + PTimesA.Z
            end
            local T3: number do --(Ax*Vx + Ay*Vy + Az*Vz)
                T3 = DeltaAcceleration:Dot(DeltaVelocity)
            end
            local T4: number do --(Ax^2 + Ay^2 + Az^2)/4
                local ASquare = DeltaAcceleration*DeltaAcceleration
                T4 = (ASquare.X+ASquare.Y+ASquare.Z)/4
            end
            if DeltaJerk == nil then
                return T0,T1,T2,T3,T4
            else
                do
                    T3 += DeltaJerk:Dot(DeltaPosition) --(Ax*Vx + Ay*Vy + Az*Vz + Jx*Px + Jy*Py + Jz*Pz)
                    T4 += DeltaJerk:Dot(DeltaVelocity)/3 --(Ax^2 + Ay^2 + Az^2)/4 + (Jx*Vx + Jy*Vy + Jz*Vz)/3
                end
                local T5: number do --(Jx*Ax + Jy*Ay + Jz*Az)/6
                    T5 = DeltaJerk:Dot(DeltaAcceleration)/6
                end
                local T6: number do --(Jx^2 + Jy^2 + Jz^2)/36
                    local JSquare = DeltaJerk*DeltaJerk
                    T6 = (JSquare.X+JSquare.Y+JSquare.Z)/36
                end
                return T0,T1,T2,T3,T4,T5,T6
            end
        end
    end
    return T0
end
--[=[
    @within BallisticsFunctions
    @private
    A workaround function to address how odd roots of a negative returns NaN when they should be returning a negative odd-rooted number.
    
    @param n -- The number of which the root will be applied to (n^root).
    @param root -- The number of the root (n^root). By default is equal to 2.
    @return number -- Returns the result of n^root. Returns NaN if the root is even and n is negative.
]=]
function BallisticUtilities:GetRoot(n: number, root: number?): number
    root = root or 2
    local NSign = math.sign(n)
    local NValue = math.abs(n)
    local ValueResult = NValue^(1/root)
    if root%2 == 1 then
        return ValueResult*NSign
    elseif NSign >= 0 then
        return ValueResult
    else
        return tonumber("nan")
    end
end
--[=[
    @within BallisticsFunctions
    @private
    A function that solves a polynomial given its coefficients up to the fourth power (T0 + T1x + T2x^2 + T3x^3 + T4x^4 = 0).
    
    @param T0 -- The first coefficient of the polynomial. T0 in T0 + T1x + T2x^2 + T3x^3 + T4x^4 = 0
    @param T1 -- The second coefficient of the polynomial. T1 in T0 + T1x + T2x^2 + T3x^3 + T4x^4 = 0
    @param T2 -- The third coefficient of the polynomial. T2 in T0 + T1x + T2x^2 + T3x^3 + T4x^4 = 0
    @param T3 -- The fourth coefficient of the polynomial. T3 in T0 + T1x + T2x^2 + T3x^3 + T4x^4 = 0
    @param T4 -- The fifth coefficient of the polynomial. T4 in T0 + T1x + T2x^2 + T3x^3 + T4x^4 = 0
    @return number, number?, number?, number? -- Returns all possible solutions of the given polynomial, up to four solutions.
]=]
function BallisticUtilities:SolvePolynomial(T0: number, T1: number, T2: number, T3: number, T4: number): (number, number|nil, number|nil, number|nil)
    T4 = T4 or 0
    T3 = T3 or 0
    T2 = T2 or 0
    T1 = T1 or 0
    T0 = T0 or 0
    local output1: number, output2: number?, output3: number?, output4: number?
    if T4 ~= 0 then --Quartic, Page 422, http://inis.jinr.ru/sl/vol1/CMC/Graphics_Gems_1,ed_A.Glassner.pdf
        --Also refer to https://mathworld.wolfram.com/QuarticEquation.html
        --A quartic equation, T4x^4 + T3x^3 + T2x^2 + T1x + T0 = 0,
        --is divided by T4: x^4 + Ax^3 + Bx^2 + Cx + D = 0
        local A,B,C,D = T3/T4,T2/T4,T1/T4,T0/T4
        --and substitution of x = y - A/4
        --eliminates the cubic term y^4 + py^2 + qy + r = 0

        --[[
            (y-A/4)^4 = y^4 - Ay^3 + 3y^2A^2/8 - xA^3/16 + A^4/256
            (y-A/4)^3 = y^3 - 3y^2A/4 + 3yA^2/16 - A^3/64
            (y-A/4)^2 = y^2 - yA/2 + A^2/16
            (y-A/4)
            0

            y^4 - Ay^3 + 3y^2A^2/8 - xA^3/16 + A^4/256 +
            A(y^3 - 3y^2A/4 + 3yA^2/16 - A^3/64) +
            B(y^2 - yA/2 + A^2/16) +
            C(y-A/4) +
            D
            y^4 - Ay^3 + 3y^2A^2/8 - xA^3/16 + A^4/256 +
            Ay^3 - 3y^2A^2/4 + 3yA^3/16 - A^4/64 +
            By^2 - yAB/2 + A^2B/16 +
            Cy-AC/4 +
            D

            y^4
            - 3y^2A^2/8 + By^2
            - yA^3/16 + 3yA^3/16 - yAB/2 + Cy
            -AC/4 + A^4/256 - A^4/64 + A^2B/16 + D

            y^4
            (-3A^2/8 + B)y^2
            (- A^3/16 + 3A^3/16 - AB/2 + C)y
            -AC/4 + A^4/256 - A^4/64 + A^2B/16 + D

            y^4
            (-3A^2/8 + B)y^2
            (2A^3/16 - AB/2 + C)y
            -AC/4 + A^4/256 - A^4/64 + A^2B/16 + D

            y^4
            (-3A^2/8 + B)y^2
            (2A^3/16 - AB/2 + C)y
            -AC/4 + 3A^4/256 + A^2B/16 + D
        ]]
        local p = -3*A*A/8 + B
        local q = A*A*A/8 - A*B/2 + C
        local r = A*C/-4 + 3*A*A*A*A/256 + A*A*B/16 + D
        --The resolvent cubic is then
        --z^3 - p/2z^2 - rz + rp/2 - q^2/8 = 0
        local z1,_z2,_z3 = self:SolvePolynomial(r*p/2-q*q/8,-1*r,p/-2,1) --in theory it should only be z1
        if z1 <= 0 then
            warn("Z is 0 or negative! This doesn't seem right...",debug.traceback())
        end
        --[[
            print(z1,_z2,_z3)
            if z2 or z3 then
                warn("Unexpected z2 and/or z3 in Quartic!")
            end
        ]]
        --With z being one root of the above equation, the roots of the quartic can
        --be obtained by solving the two quadratic equations
        --y^2 +- y*sqrt(2z-p) + z -+ sqrt(z^2-r) = 0
        
        --[[
            y^2 + y*sqrt(2z-p) + z - sqrt(z^2-r) = 0
            y^2 - y*sqrt(2z-p) + z + sqrt(z^2-r) = 0
        ]]
        local SqrtZ2R = math.sqrt(z1*z1-r)
        local Sqrt2ZP = math.sqrt(2*z1-p)
        local y1,y2 = self:SolvePolynomial(z1 - SqrtZ2R, Sqrt2ZP, 1)
        local y3,y4 = self:SolvePolynomial(z1 + SqrtZ2R, -1*Sqrt2ZP, 1)
        --Resubstitution yields the correct values for x.
        if y1 ~= nil then
            output1 = y1 - A/4
        end
        if y2 ~= nil then
            output2 = y2 - A/4
        end
        if y3 ~= nil then
            output3 = y3 - A/4
        end
        if y4 ~= nil then
            output4 = y4 - A/4
        end
        --output1,output2,output3,output4 = y1,y2,y3,y4
    elseif T3 ~= 0 then --Cubic, Page 420, http://inis.jinr.ru/sl/vol1/CMC/Graphics_Gems_1,ed_A.Glassner.pdf
        --A quartic equation, T3x^3 + T2x^2 + T1x + T0 = 0,
        --is divided by T3: x^3 + Ax^2 + Bx + C = 0
        local A,B,C = T2/T3,T1/T3,T0/T3
        --print("A:",A)
        --and substitution of x = y - A/3
        --eliminates the cubic term y^3 + 3py + 2q = 0

        --[[
            (y - A/3)^3 = y^3 - A*y^2 + A^2*y/3 - A^3/27
            (y - A/3)^2 = y^2 - 2*Ay/3 + A^2/9
            (y - A/3) = (y - A/3)
            0
            
            y^3 - A*y^2 + A^2*y/3 - A^3/27
            + A*y^2 - 2*A^2*y/3 + A^3/9
            + B*y - A*B/3
            + C

            y^3
            - A*y^2 + A*y^2
            - 2*A^2*y/3 + A^2*y/3 + B*y
            - A^3/27 + A^3/9 - A*B/3 + C

            y^3
            (-*A^2/3 + A^2/3 + B)*y
            - A^3/27 + A^3/9 - A*B/3 + C

            y^3
            (B - *A^2/3)*y
            - A^3/27 + 3*A^3/27 - A*B/3 + C

            y^3
            (B - *A^2/3)*y
            2*A^3/27 - A*B/3 + C

            y^3
            (B/3 - A^2/9)*y*3
            (A^3/27 - A*B/6 + C/2)*2
        ]]
        local p = B/3 - A*A/9
        local q = A*A*A/27 - A*B/6 + C/2
        --Using Cardano's Formula (G. Cardano, 1501,1576), the determinant is
        -- D = q^2 + p^3
        --u,v = (-q +- sqrt(D))^(1/3)
        --and the roots are
        --y1 = u+v
        --y2,3 = -(u+v)/2 +- sqrt(3)/2*(u-v)i
        local PCubed = p*p*p
        local D = q*q + PCubed
        --[[Three cases can be distinguished:
            D > 0: one real (yl
            ), two conjugated complex values (y3, y3)
            D = 0: two real values, y2 = y3
            D < 0: three different real values.
        ]]
        local y1:number?,y2:number?,y3:number?
        if D >= 0 then --At least one root, potentially two
            local SqrtD = math.sqrt(D)
            local u = self:GetRoot(-1*q + SqrtD,3)
            --print("u:",(-1*q + SqrtD),u)
            local v = self:GetRoot(-1*q - SqrtD,3)
            --print("v:",(-1*q - SqrtD),v)
            y1 = u+v
            --print("y1:",y1)
            if D < 1e-6 then --1e-6 is essentially 0 due to float errors, two roots
                y2 = -(y1)/2
                --print("y2:",y2)
            end
        else --Three roots
            --[[
                Messy solution:
                u = (-q + math.sqrt(D))^(1/3)
                v = (-q - math.sqrt(D))^(1/3)
                output2 = -(u+v)/2 + math.sqrt(3)/2*(u-v)i
                if D < 0 then
                    output3 = -(u+v)/2 - math.sqrt(3)/2*(u-v)i
                end
            ]]
            --[[
                In the case of D < 0 (the so-called casus irreducibilis) trigonometric
                substitution helps to find all three solutions without the need for complex
                arithmetics:
            ]]
            -- cosR = -q/sqrt(-p^3)
            -- y1 = 2*sqrt(-p)*cos(R/3)
            -- y2,3 = -2*sqrt(-p)*cos((R+-pi)/3)
            
            local R = math.acos(-1*q/math.sqrt(-1*PCubed))
            local NSqrtP = math.sqrt(-1*p)
            y1 = 2*NSqrtP*math.cos(R/3)
            y2 = -2*NSqrtP*math.cos((R+math.pi)/3)
            y3 = -2*NSqrtP*math.cos((R-math.pi)/3)
        end
        --Resubstitution yields the correct values for x.
        if y1 ~= nil and y1 == y1 then
            output1 = y1 - A/3
        end
        if y2 ~= nil and y2 == y2 then
            output2 = y2 - A/3
        end
        if y3 ~= nil and y3 == y3 then
            output3 = y3 - A/3
        end
        --output1,output2,output3 = y1,y2,y3
    elseif T2 ~= 0 then -- Quadratic Returns how many solutions there are, https://www.forrestthewoods.com/blog/solving_ballistic_trajectories/
        -- 0 = T0 + T1v + T2v^2
        if math.abs(T2)<1e-6 then
            if math.abs(T1)<1e-6 then
                if math.abs(T0)<1e-6 then
                    output1 = 0
                end
            else
                output1 = -T0/T1
            end
        else
            local D = T1^2-4*T2*T0
            local sqrtD = math.sqrt(D)
            if D >= 0 then
                output1 = (-T1 + sqrtD) / (2*T2)
                if D > 0 then
                    output2 = (-T1 - sqrtD) / (2*T2)
                end
            end
        end
    elseif T1 ~= 0 then
        output1 = -T0/T1
    else
        output1 = T0
    end
    --print("SolvePolynomial Inputs:",T0,T1,T2,T3,T4)
    --print("SolvePolynomial Outputs:",output1,output2,output3,output4)
    return output1,output2,output3,output4
end
--[=[
    @within BallisticsFunctions
    @private
    A short function that finds the overall vector between two Vector3s (Final Vector - Initial Vector).
    
    @param Target -- The vector3 of the target. By default Vector3.new(0,0,0).
    @param Shooter -- The vector3 of the shooter. By default Vector3.new(0,0,0).
    @return Vector3 -- The resulting Vector3 of (Target - Shooter).
]=]
function BallisticUtilities:ProduceDeltas(Target: Vector3?, Shooter: Vector3?): Vector3
    if Target or Shooter then
        Target = Target or Vector3.new()
        Shooter = Shooter or Vector3.new()
        return Target - Shooter
    end
end
--[=[
    @within BallisticsFunctions
    @private
    A function that returns the derivative of a polynomial, given its coefficients.
    
    @param ... -- The coefficients of the polynomial. T0, T1, T2, T3, T4, ..., Tn, where: T0 + T1x + T2x^2 + T3x^3 + T4x^4 + ... + Tn^n = 0
    @return ...number -- The coefficients of the derviative polynomial from the initial polynomial.
]=]
function BallisticUtilities:ProduceDerivative(...:number): ...number
    local Coefficients: Array<number> = table.pack(...)
    local NewCoefficients = {}
    for i=2, #Coefficients do --skip 1 as that'll just equal 0
        local ThisCoefficient = Coefficients[i]
        local ThisPower = i-1
        table.insert(NewCoefficients,ThisCoefficient*ThisPower)
    end
    return table.unpack(NewCoefficients)
end
--[=[
    @within BallisticsFunctions
    @private
    A function that returns the result of pluging in Input into a polynomial, given its coefficients. It can also evaluate the limit of positive and negative math.huge (infinity) if the Input is set to be that.
    
    @param Input -- The number that will be put into the given polynomial. Can be math.huge (infinity).
    @param ... -- The coefficients of the polynomial. T0, T1, T2, T3, T4, ..., Tn, where: T0 + T1x + T2x^2 + T3x^3 + T4x^4 + ... + Tn^n = 0
    @return number -- The result of plugging in Input into the given polynomial.
]=]
function BallisticUtilities:InputPolynomial(Input: number, ...:number): number
    local Coefficients: Array<number> = table.pack(...)
    if math.abs(Input) == math.huge then
        --task.wait(0.1)
        local Sign = math.sign(Input) -- -1
        local Power = #Coefficients-1 -- 5
        local LargestCoefficient = Coefficients[Power+1]
        if Power%2 == 0 then
            Sign = 1
        end
        Sign *= math.sign(LargestCoefficient)
        return Sign*math.huge -- math.sign(InputPolynomial(1,...))*math.huge --
    end
    local Sum = 0
    for i=1, #Coefficients do
        local ThisCoefficient = Coefficients[i]
        local ThisPower = i-1
        Sum += ThisCoefficient*(Input^(ThisPower))
    end
    return Sum
end
--[=[
    @within BallisticsFunctions
    @private
    A function that returnsif two numbers are within a given range of precision.
    
    @param N1 -- The first number.
    @param N2 -- The second number.
    @param Precision -- The desired range of precision that the numbers should be in. This number will be divided by 2 and added/subtracted to the first number to check to see if the second number is inbetween it. By default is set to [BallisticUtilities.Precision].
    @return boolean -- The result of plugging in Input into the given polynomial.
]=]
function BallisticUtilities:CompareNumbers(N1: number, N2: number, Precision: number?): boolean
    Precision = Precision or self.Precision
    if N1 == N2 then
        return true
    end
    local N1L = N1 - Precision/2
    local N1G = N1 + Precision/2
    return (N1L <= N2 and N2 <= N1G)
end
--[=[
    @within BallisticsFunctions
    @private
    A function that recursively calls itself to figure out where inbetween Point1 and Point2 it equals 0 within the range of [BallisticUtilities.Precision]. It performs the Bisection method to do this. If either point is equal to math.huge (infinity), it will replace it with [BallisticUtilities.MaxNumber].
    
    @param Point1 -- One critical point of which should be a min/max.
    @param Point2 -- Another critical poit of which should be a min/max.
    @param ... -- The coefficients of the polynomial. T0, T1, T2, T3, T4, ..., Tn, where: T0 + T1x + T2x^2 + T3x^3 + T4x^4 + ... + Tn^n = 0
    @return number -- The result of plugging in Input into the given polynomial.
]=]
function BallisticUtilities:ProduceEstimate(Point1: number, Point2: number, ...:number): number
    --print("Produce Estimate Between", Point1, "and", Point2, "for",...)
    local Midpoint: number
    local Point1Estimate,Point2Estimate = self:InputPolynomial(Point1,...),self:InputPolynomial(Point2,...)
    local LowerBound:number, HigherBound:number, LowerPoint:number, HigherPoint:number do
        if Point1Estimate < Point2Estimate then
            --print("Point1Estimate < Point2Estimate")
            LowerBound = Point1Estimate
            LowerPoint = Point1
            HigherBound = Point2Estimate
            HigherPoint = Point2
        else
            --print("Point1Estimate > Point2Estimate")
            LowerBound = Point2Estimate
            LowerPoint = Point2
            HigherBound = Point1Estimate
            HigherPoint = Point1
        end
    end
    --print(LowerBound,"<",HigherBound)
    local LowerAbs, HigherAbs = math.abs(LowerBound), math.abs(HigherBound)
    --print(LowerAbs,"<",HigherAbs)
    if LowerAbs == math.huge or HigherAbs == math.huge then
        --print("Inf trigger")
        if LowerAbs ~= HigherAbs then
            if LowerAbs == math.huge then
                --print("LowerAbs == math.huge")
                LowerBound = -1 * self.MaxNumber
                Midpoint = math.sign(LowerPoint) * self.MaxNumber
            else
                --print("HigherAbs == math.huge")
                HigherBound = self.MaxNumber
                Midpoint = math.sign(HigherPoint) * self.MaxNumber
            end
        else
            Midpoint = 0
        end
    else
        Midpoint = (Point1 + Point2)/2
    end
    --print("Midpoint:",Midpoint)
    local CurrentEstimate = self:InputPolynomial(Midpoint,...)
    --print("Current Estimate Ranges:",LowerBound,"-",CurrentEstimate,"-",HigherBound)
    if LowerBound > 0 and HigherBound < 0 then
        warn("No solutions that equal 0 between points",Point1,"and",Point2,"for coefficients",...,".")
        return nil
    elseif self:CompareNumbers(0, CurrentEstimate) then
        --print("Midpoint Generated [New Point, Result For Error]:",Midpoint,",",CurrentEstimate)
        return Midpoint
    else
        if LowerBound < CurrentEstimate and CurrentEstimate < HigherBound then
            if CurrentEstimate < 0 then
                return self:ProduceEstimate(Midpoint, HigherPoint, ...)
            else
                return self:ProduceEstimate(LowerPoint, Midpoint, ...)
            end
        elseif LowerBound < CurrentEstimate then
            return self:ProduceEstimate(LowerPoint, Midpoint, ...)
        elseif CurrentEstimate < HigherBound then
            return self:ProduceEstimate(Midpoint, HigherPoint, ...)
        else
            warn("Cannot find valid midpoint between points",Point1,"and",Point2,"for coefficients",...,"!")
            return nil
        end
    end
end
--[=[
    @within BallisticsFunctions
    @private
    A function that creates a new array which only contains unique values.
    
    @param InputArray -- The array of which should be pruned of duplicates.
    @return Array<any> -- The result of plugging in Input into the given polynomial.
]=]
function BallisticUtilities:RemoveDuplicatesFromArray(InputArray: Array<any>): Array<any>
    local output = {}
    for i=1, #InputArray do
        local v = InputArray[i]
        if not table.find(output,v) then
            table.insert(output,v)
        end
    end
    return output
end
--[=[
    @within BallisticsFunctions
    @private
    A function that finds and returns the values where a relative/absolute max/min occurs in the polynomial.
    If the polynomial is a quartic or lower, it will use the general formula in SolvePolynomial.
    Otherwise, it will perform the bisection method via ProduceEstimate.
    
    @param ... -- The coefficients of the polynomial. T0, T1, T2, T3, T4, ..., Tn, where: T0 + T1x + T2x^2 + T3x^3 + T4x^4 + ... + Tn^n = 0
    @return ...number -- The critical points of the given polynomial.
]=]
function BallisticUtilities:ProduceCriticalPoints(...:number): ...number
    print("Produce Critical Points for",...)
    local Coefficients: Array<number> = table.pack(...)
    local HighestPower = #Coefficients-1
    local CriticalPoints: Array<number>
    if HighestPower > 4 then --needs a check to see if enough solutions are given, if not, go deeper for derivatives
        local DerivativeCriticalPoints = table.pack(self:ProduceCriticalPoints(self:ProduceDerivative(...)))
        DerivativeCriticalPoints = self:RemoveDuplicatesFromArray(DerivativeCriticalPoints)
        table.sort(DerivativeCriticalPoints,function(a,b)
            return a < b
        end)
        --print(table.unpack(DerivativeCriticalPoints))
        CriticalPoints = {-math.huge,math.huge}
        for i=2, #DerivativeCriticalPoints do
            local PreviousPoint = DerivativeCriticalPoints[i-1]
            local ThisPoint = DerivativeCriticalPoints[i]
            local PotentialPoint = self:ProduceEstimate(PreviousPoint,ThisPoint,...)
            if PotentialPoint ~= nil then
                table.insert(CriticalPoints,PotentialPoint)
            end
        end
    else
        CriticalPoints = table.pack(-math.huge,math.huge,self:SolvePolynomial(...))
    end
    table.sort(CriticalPoints, function(a:number,b:number)
        return a < b
    end)
    print("Critical Points for",...,":",table.unpack(CriticalPoints))
    return table.unpack(CriticalPoints,1,#CriticalPoints)
end
--[=[
    @within BallisticsFunctions
    @private
    A function that finds and returns the solutions to the inputted polynomial which doesn't have a general formula with critical points and the bisection method.
    
    @param ... -- The coefficients of the polynomial that needs to be solved. T0, T1, T2, T3, T4, ..., Tn, where: T0 + T1x + T2x^2 + T3x^3 + T4x^4 + ... + Tn^n = 0
    @return ...number -- The solutions of the given polynomial.
]=]
function BallisticUtilities:GetEstimate(...:number): ...number
    local PotentialSolutions = table.pack(self:ProduceCriticalPoints(...))
    local ActualSolutions = {}
    for i=1, #PotentialSolutions do
        local ThisSolution = PotentialSolutions[i]
        --print(ThisSolution)
        if self:CompareNumbers(self:InputPolynomial(ThisSolution,...),0) then
            table.insert(ActualSolutions, ThisSolution)
        end
    end
    return table.unpack(ActualSolutions,1,#ActualSolutions)
end
--[=[
    @class BallisticsFunctions

    The public table that contains the functions for finding the points of collision between two projectiles.
]=]
--[=[
    @prop Utilities BallisticUtilities
    @private
    @within BallisticsFunctions
    The BallisticUtilities class of which contains all of the nessessary functions used to create the right data.
]=]
local BallisticsFunctions = {
    Utilities = BallisticUtilities;
} do
    BallisticsFunctions.__index = BallisticsFunctions
    BallisticsFunctions = setmetatable(BallisticsFunctions,BallisticsFunctions)
end
--[=[
    @within BallisticsFunctions
    @type LookVector Vector3
    A unit Vector3 that is equivalent to CFrame.LookVector.
]=]
--[=[
    @within BallisticsFunctions
    @type AverageHitPosition Vector3
    A Vector3 where the projectile and the target will collide.
]=]
--[=[
    @within BallisticsFunctions
    @type PositionAccuracy number
    A number that shows how far off the average is from the projected target collision positions due to calculation incaccuracies.
]=]
--[=[
    @within BallisticsFunctions
    @type SimulatedHitPosition Vector3
    The projected target collision position via: HitPosition = TargetInitialPosition + TargetVelocity*TravelTime + TargetAcceleration/2*TravelTime^2 + TargetJerk/3*TravelTime^3.
]=]
--[=[
    @within BallisticsFunctions
    @type ProjectedHitPosition Vector3
    The projected projectile collision position via: HitPosition = ShooterInitialPosition + ((SimulationLookVector*ProjectileSpeed + ShooterVelocity)*TravelTime) + ShooterAcceleration/2*TravelTime^2 + ShooterJerk/3*TravelTime^3.
]=]
--[=[
    @private
    @within BallisticsFunctions
    An internal function that computes the LookVector and HitPosition from the given results and variables.
    
    @param results -- An array that contains solutions from either GetHitTimes or GetHitTimesWithJerk.
    @param ProjectileSpeed -- The initial speed of the projectile.
    @param ShooterPosition -- The (initial) position of the projectile/shooter.
    @param ShooterVelocity -- The (initial) velocity of the projectile/shooter.
    @param ShooterAcceleration -- The (initial) acceleration of the projectile/shooter.
    @param ShooterJerk -- The jerk (rate of change of acceleration) of the projectile/shooter.
    @param TargetPosition -- The (initial) position of the target.
    @param TargetVelocity -- The (initial) velocity of the target.
    @param TargetAcceleration -- The (initial) velocity of the target.
    @param TargetJerk -- The jerk (rate of change of acceleration) of the target.
    @return LookVector?, AverageHitPosition?, PositionAccuracy?, SimulatedHitPosition?, ProjectedHitPosition? -- Returns nil if there isn't a valid hit direction.
]=]
local function ReturnHitInfo(results: Array<number>, ProjectileSpeed: number, ShooterPosition: Vector3, ShooterVelocity: Vector3?, ShooterAcceleration: Vector3?, ShooterJerk: Vector3?, TargetPosition: Vector3, TargetVelocity: Vector3?, TargetAcceleration: Vector3?, TargetJerk: Vector3?): (Vector3, Vector3, number, Vector3, Vector3)
    if #results > 0 then
        --print("Hitable Check Results:", table.unpack(results))
        table.sort(results,function(a,b)
            if a >= 0 then
                if b >= 0 then
                    return a < b
                else
                    return true
                end
            else
                return false
            end
        end)
        local MinimalTime = results[1]
        --print("Minimal Time to hit the target:",MinimalTime)
        assert(MinimalTime ~= nil and MinimalTime >= 0, "Minimal Time is negative or doesn't exist!")
        local ProjectedHitPosition = TargetPosition + TargetVelocity * MinimalTime + TargetAcceleration * MinimalTime*MinimalTime / 2 + TargetJerk * MinimalTime*MinimalTime*MinimalTime / 3
        --print("Target Intercept Position:", ProjectedHitPosition)
        local SimulationLookVector = (ProjectedHitPosition - ShooterPosition).Unit
        --print("Simulated Look Vector to hit target:", SimulationLookVector)
        local SimulatedNetVelocity = (SimulationLookVector * ProjectileSpeed + ShooterVelocity) * MinimalTime
        local SimulatedHitPosition = ShooterPosition + SimulatedNetVelocity + ShooterAcceleration * MinimalTime*MinimalTime / 2 + ShooterJerk * MinimalTime*MinimalTime*MinimalTime / 3
        --print("Simulated Hit Position:", SimulatedHitPosition)
        local AverageHitPosition = (SimulatedHitPosition + ProjectedHitPosition)/2
        local LargestMagnitudeVector:Vector3 do
            local ProjectedMagnitude = ProjectedHitPosition.Magnitude
            local SimulatedMagnitude = SimulatedHitPosition.Magnitude
            if ProjectedMagnitude > SimulatedMagnitude then
                LargestMagnitudeVector = ProjectedHitPosition
            else
                LargestMagnitudeVector = SimulatedHitPosition
            end
        end
        local PositionAccuracy = LargestMagnitudeVector-AverageHitPosition
        return SimulationLookVector, AverageHitPosition, PositionAccuracy, SimulatedHitPosition, ProjectedHitPosition
    else
        --warn("Hitable Check failure! Not enough numbers were returned!")
        return nil
    end
end
--[=[
    @since v0.2.1
    @tag Untested
    A function that computes the LookVector and HitPosition from ProjectileSpeed, Positions, Velocities, and Acceleration. Refer to [BallisticsFunctions:GetHitInfoWithJerk] to compute with Jerk.
    
    @param ProjectileSpeed -- The initial speed of the projectile.
    @param ShooterPosition -- The (initial) position of the projectile/shooter.
    @param ShooterVelocity -- The (initial) velocity of the projectile/shooter.
    @param ShooterAcceleration -- The (initial) acceleration of the projectile/shooter.
    @param TargetPosition -- The (initial) position of the target.
    @param TargetVelocity -- The (initial) velocity of the target.
    @param TargetAcceleration -- The (initial) velocity of the target.
    @return LookVector?, AverageHitPosition?, PositionAccuracy?, SimulatedHitPosition?, ProjectedHitPosition? -- Returns nil if there isn't a valid hit direction.
]=]
function BallisticsFunctions:GetHitInfo(ProjectileSpeed: number, ShooterPosition: Vector3, ShooterVelocity: Vector3?, ShooterAcceleration: Vector3?, TargetPosition: Vector3, TargetVelocity: Vector3?, TargetAcceleration: Vector3?): (Vector3, Vector3, number, Vector3, Vector3)
    local results = table.pack(self:GetHitTimes(ProjectileSpeed, ShooterPosition, ShooterVelocity, ShooterAcceleration, TargetPosition, TargetVelocity, TargetAcceleration))
    return ReturnHitInfo(results, ProjectileSpeed, ShooterPosition, ShooterVelocity, ShooterAcceleration, nil, TargetPosition, TargetVelocity, TargetAcceleration, nil)
end
--[=[
    @since v0.2.1
    @tag Untested
    A function that computes the LookVector and HitPosition from ProjectileSpeed, Positions, Velocities, Acceleration, and Jerk.
    
    @param ProjectileSpeed -- The initial speed of the projectile.
    @param ShooterPosition -- The (initial) position of the projectile/shooter.
    @param ShooterVelocity -- The (initial) velocity of the projectile/shooter.
    @param ShooterAcceleration -- The (initial) acceleration of the projectile/shooter.
    @param ShooterJerk -- The jerk (rate of change of acceleration) of the projectile/shooter.
    @param TargetPosition -- The (initial) position of the target.
    @param TargetVelocity -- The (initial) velocity of the target.
    @param TargetAcceleration -- The (initial) velocity of the target.
    @param TargetJerk -- The jerk (rate of change of acceleration) of the target.
    @return LookVector?, AverageHitPosition?, PositionAccuracy?, SimulatedHitPosition?, ProjectedHitPosition? -- Returns nil if there isn't a valid hit direction.
]=]
function BallisticsFunctions:GetHitInfoWithJerk(ProjectileSpeed: number, ShooterPosition: Vector3, ShooterVelocity: Vector3?, ShooterAcceleration: Vector3?, ShooterJerk: Vector3?, TargetPosition: Vector3, TargetVelocity: Vector3?, TargetAcceleration: Vector3?, TargetJerk: Vector3?): (Vector3, Vector3, number, Vector3, Vector3)
    local results = table.pack(self:GetHitTimesWithJerk(ProjectileSpeed, ShooterPosition, ShooterVelocity, ShooterAcceleration, ShooterJerk, TargetPosition, TargetVelocity, TargetAcceleration, TargetJerk))
    return ReturnHitInfo(results, ProjectileSpeed, ShooterPosition, ShooterVelocity, ShooterAcceleration, ShooterJerk, TargetPosition, TargetVelocity, TargetAcceleration, TargetJerk)
end

--[[
function BallisticsFunctions:GetHitTimesWithDictionary(Input: Dictionary<any>)
    local ProjectileSpeed: number? = Input.ProjectileSpeed
    assert(type(ProjectileSpeed) == "number", "Input.ProjectileSpeed is not a number! Debug:"..debug.traceback())
    local ConserveMomentum: boolean? = Input.ConserveMomentum or false
    local ShooterPosition: Vector3? = Input.ShooterPosition
    assert(typeof(ShooterPosition) == "Vector3", "Input.ShooterPosition is not a Vector3! Debug:"..debug.traceback())
    local ShooterVelocity: Vector3?,ShooterAcceleration: Vector3?
    if ConserveMomentum then
        ShooterVelocity = Input.ShooterVelocity
        ShooterAcceleration = Input.ShooterAcceleration
    else
        ShooterVelocity = Vector3.new()
        ShooterAcceleration = Vector3.new()
    end
    local TargetPosition: Vector3? = Input.TargetPosition
    local TargetVelocity: Vector3? = Input.TargetVelocity
    local TargetAcceleration: Vector3? = Input.TargetAcceleration
end
]]

--[=[
    @since v0.2.0
    A function that computes the times of collision from ProjectileSpeed, Positions, Velocities, and Acceleration. Refer to [BallisticsFunctions:GetHitTimesWithJerk] to compute with Jerk.
    
    @param ProjectileSpeed -- The initial speed of the projectile.
    @param ShooterPosition -- The (initial) position of the projectile/shooter.
    @param ShooterVelocity -- The (initial) velocity of the projectile/shooter.
    @param ShooterAcceleration -- The (initial) acceleration of the projectile/shooter.
    @param TargetPosition -- The (initial) position of the target.
    @param TargetVelocity -- The (initial) velocity of the target.
    @param TargetAcceleration -- The (initial) velocity of the target.
    @return ...number -- The times of when the projectile and target can collide.
]=]
function BallisticsFunctions:GetHitTimes(ProjectileSpeed: number, ShooterPosition: Vector3, ShooterVelocity: Vector3?, ShooterAcceleration: Vector3?, TargetPosition: Vector3, TargetVelocity: Vector3?, TargetAcceleration: Vector3?)
    local DeltaPosition = self.Utilities:ProduceDeltas(TargetPosition, ShooterPosition) --p
    local DeltaVelocity = self.Utilities:ProduceDeltas(TargetVelocity, ShooterVelocity) --v
    local DeltaAcceleration = self.Utilities:ProduceDeltas(TargetAcceleration, ShooterAcceleration) do --a
        if DeltaAcceleration then
            DeltaVelocity = DeltaVelocity or Vector3.new()
        end
    end
    local output = table.pack(self.Utilities:SolvePolynomial(self.Utilities:ProduceCoefficients(ProjectileSpeed, DeltaPosition, DeltaVelocity, DeltaAcceleration)))
    --[[for i=1, #output do
        local OriginalOutput = output[i]
        output[i] = math.round(OriginalOutput/self.Utilities.Precision)*self.Utilities.Precision
    end]]
    --print("GetHitTimes Output:",table.unpack(output))
    return table.unpack(output)
end
--[=[
    @since v0.2.0
    A function that computes the times of collision from ProjectileSpeed, Positions, Velocities, Acceleration, and Jerk.
    
    @param ProjectileSpeed -- The initial speed of the projectile.
    @param ShooterPosition -- The (initial) position of the projectile/shooter.
    @param ShooterVelocity -- The (initial) velocity of the projectile/shooter.
    @param ShooterAcceleration -- The (initial) acceleration of the projectile/shooter.
    @param ShooterJerk -- The jerk (rate of change of acceleration) of the projectile/shooter.
    @param TargetPosition -- The (initial) position of the target.
    @param TargetVelocity -- The (initial) velocity of the target.
    @param TargetAcceleration -- The (initial) velocity of the target.
    @param TargetJerk -- The jerk (rate of change of acceleration) of the target.
    @return ...number -- The times of when the projectile and target can collide.
]=]
function BallisticsFunctions:GetHitTimesWithJerk(ProjectileSpeed: number, ShooterPosition: Vector3, ShooterVelocity: Vector3?, ShooterAcceleration: Vector3?, ShooterJerk: Vector3?, TargetPosition: Vector3, TargetVelocity: Vector3?, TargetAcceleration: Vector3?, TargetJerk: Vector3?)
    local DeltaPosition = self.Utilities:ProduceDeltas(TargetPosition, ShooterPosition) --p
    local DeltaVelocity = self.Utilities:ProduceDeltas(TargetVelocity, ShooterVelocity) --v
    local DeltaAcceleration = self.Utilities:ProduceDeltas(TargetAcceleration, ShooterAcceleration) do --a
        if DeltaAcceleration then
            DeltaVelocity = DeltaVelocity or Vector3.new()
        end
    end
    local DeltaJerk = self.Utilities:ProduceDeltas(TargetJerk, ShooterJerk) do --j
        if DeltaJerk then
            DeltaVelocity = DeltaVelocity or Vector3.new()
            DeltaAcceleration = DeltaAcceleration or Vector3.new()
        end
    end
    local output do
        if DeltaJerk ~= nil then
            --print(self.Utilities:ProduceCoefficients(ProjectileSpeed, DeltaPosition, DeltaVelocity, DeltaAcceleration, DeltaJerk))
            output = table.pack(self.Utilities:GetEstimate(self.Utilities:ProduceCoefficients(ProjectileSpeed, DeltaPosition, DeltaVelocity, DeltaAcceleration, DeltaJerk)))
        else
            output = table.pack(self.Utilities:SolvePolynomial(self.Utilities:ProduceCoefficients(ProjectileSpeed, DeltaPosition, DeltaVelocity, DeltaAcceleration)))
        end
    end
    --[[for i=1, #output do
        local OriginalOutput = output[i]
        output[i] = math.round(OriginalOutput/self.Utilities.Precision)*self.Utilities.Precision
    end]]
    --print("GetHitTimes Output:",table.unpack(output))
    return table.unpack(output)
end
return BallisticsFunctions