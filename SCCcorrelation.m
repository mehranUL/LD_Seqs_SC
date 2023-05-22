function SCC = SCCcorrelation(Control1, Control2, N)
               
                a = 0;
                b = 0;
                c = 0;
                d = 0;

                for z = 1:N
                    if Control1(z) == 1 && Control2(z) == 1
                        a = a + 1;
                    end
                    if Control1(z) == 1 && Control2(z) == 0
                        b = b + 1;
                    end
                    if Control1(z) == 0 && Control2(z) == 1
                        c = c + 1;
                    end
                    if Control1(z) == 0 && Control2(z) == 0
                        d = d + 1;
                    end
                end

                if (a*d) > (b*c)
                    SCC = ((a*d)-(b*c))/((N*(min((a+b),(a+c))))-((a+b)*(a+c)));
                end
                if (a*d) <= (b*c)
                    SCC = ((a*d)-(b*c))/(((a+b)*(a+c))-(N*(max((a-d),0))));
                end

                if isnan(SCC)
                    SCC = 0;
                end
end