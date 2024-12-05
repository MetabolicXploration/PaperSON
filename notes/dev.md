## v1 Design ideas

### raw

'raw' it is a free structure level. Its only goal is the introduction of the data as coherent with the paper as posible. That means that by reading the paper you should be able to underestand what the 'raw.json' file contain. This includes the ids of variables, thee units, the source, etc. 

For instance:
```
{
    "table2" : {
        "F1" : {
            "vals" : [...],
            "unit" : "hours"
        }
    }
}
```

should contain the date of 'table2' of the manuscript, 'F1' should be the same symbol used in the table, and "hours" should be the units expressed in the table as well.  

#### About nesting depth

It is recommended to use nesting only when it is necessary for making the data structure clearer. Also, try not to mix different 'kind' of data, nesting can be used to separate them. Anyway, remember, 'raw' data is structure free. 

For instance:

```
{
    "Supp2" : {
        "meta" : {
            "source" : "Sumpplementary file 00118s2.xlsx"
        }, 
        "data" : {
            "D" : {
                "unit" : "1/h",
                "val" : [0.1, 0.2, 0.3, 0.4],
                "err" : [0.0, 0.0, 0.0, 0.0]
            },
        }
    }
}
```