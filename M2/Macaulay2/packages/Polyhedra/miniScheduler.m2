protect targetProperty
protect sourceProperties

rules = new CacheTable


findOneStepComputableProperties = method()
findOneStepComputableProperties Set := givenProperties -> (
   allRules := keys rules;
   goodRules := select(allRules, rule -> rule#0 <= givenProperties);
   result := apply(goodRules, rule -> rule#1);
   result = sum result;
   if result === 0 then return set {} else result
)


findRulesForProperty = method()
findRulesForProperty(Set, Set, Symbol) := (given, forbidden, target) -> (
   allRules := keys rules;
   goodRules := select(allRules, rule -> member(target, rule#1));
   allowedRules := select(goodRules, rule -> #(rule#0 * forbidden) == 0);
   return allowedRules
)

findOneStepRulesForProperty = method()
findOneStepRulesForProperty(Set, Set, Symbol) := (given, forbidden, target) -> (
   allowedRules := findRulesForProperty(given, forbidden, target);
   satisfiableRules := select(allowedRules, rule -> rule#0 <= given);
   return satisfiableRules
)

resolveTarget = method()
resolveTarget(Set, Set, Symbol) := (given, forbidden, target) -> (
   allowedRules := findRulesForProperty(given, forbidden, target);
   for rule in allowedRules do (
      resolvent := rule#0;
      newForbidden := forbidden + set {target};
      << "Resolvent: " << resolvent << endl;
      result := findPathRecursively(given, newForbidden, resolvent);
      if #result > 0 then return append(result, (rule))
   );
   return ()
)


findPathRecursively = method()
findPathRecursively(Set, Set, Set) := (given, forbidden, targets) -> (
   computableProperties := findOneStepComputableProperties(given);
   << "Computable: " << endl;
   goodTargets := targets * computableProperties;
   result := toSequence apply(toList goodTargets, target -> (findOneStepRulesForProperty(given, forbidden, target))#0);
   newGiven := given + goodTargets;
   badTargets := targets - newGiven;
   for target in toList badTargets do (
      newRules := resolveTarget(newGiven, forbidden, target);
      if #newRules == 0 then return () 
      else (
         result = result | newRules;
         newGiven = newGiven + set {target};
      )
   );
   return result
)



findSuitableRulePath = method()
findSuitableRulePath(HashTable, Set) := (HT, targets) -> (
   givenProperties := set flatten {keys HT, keys HT.cache};
   path := findPathRecursively(givenProperties, set {}, targets);
   return path
)

applySingleRule = method()
applySingleRule(HashTable, Sequence, Set) := (HT, rule, given) -> (
   if rule#1 <= given then return given
   else (
      << "Applying Rule:" << endl;
      << "Source: " << rule#0 << endl;
      << "Target: " << rule#1 << endl;
      rules#rule HT;
      given + rule#1
   )
)

applySuitableRules = method()
applySuitableRules(HashTable, Set) := (HT, targets) -> (
   givenProperties := set flatten {keys HT, keys HT.cache};
   rulePath := findSuitableRulePath(HT, targets);
   if #rulePath == 0 then error "Cannot compute properties."
   else (
      for edge in rulePath do (
         givenProperties = applySingleRule(HT, edge, givenProperties);
      )
   )
)

applySuitableRules(HashTable, Symbol) := (HT, target) -> (
   applySuitableRules(HT, set {target});
   HT.cache#target
)

makeRule = method()
makeRule(Set, Set, FunctionClosure) := (source, target, f) -> (
   key := (source, target);
   rules#key = f;
)
