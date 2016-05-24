raw = value Core#"private dictionary"#"raw"
rawPolymakeConvexHull = value Core#"private dictionary"#"rawPolymakeConvexHull"

export {
   "polymakeConvexHull"
}

polymakeConvexHull = method()
polymakeConvexHull MutableMatrix := V -> (
   rawPolymakeConvexHull(raw V)
)
