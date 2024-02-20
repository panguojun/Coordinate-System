using UnityEngine;
using System;

public struct Life
{
    public int lifespan; // 寿命
    public string status; // 状态（活着/死亡）
    public long boundaryStart; // 有效时间段的起始时间（毫秒）
    public long boundaryEnd; // 有效时间段的结束时间（毫秒）

    public Life(int lifespan, long boundaryStart, long boundaryEnd)
    {
        this.lifespan = lifespan;
        this.boundaryStart = boundaryStart;
        this.boundaryEnd = boundaryEnd;
        this.status = "alive";
    }

    public void UpdateStatus(long currentTime)
    {
         if (currentTime < boundaryStart || currentTime > boundaryEnd)
        {
            this.status = "dead";
            lifespan ++;
        }
        else
        {
            this.status = "alive";
        }
    }

    public string GetStatus()
    {
        return this.status;
    }

    // 重载加法运算符
    public static Life operator +(Life life, bool value)
    {
        if (value)
        {
            life.boundaryStart = life.boundaryStart+(1);
            life.boundaryEnd = life.boundaryEnd+(1);
        }
        else
        {
            life.boundaryStart = life.boundaryStart+(-1);
            life.boundaryEnd = life.boundaryEnd+(-1);
        }

        return life;
    }

    // 重载减法运算符
    public static Life operator -(Life life, bool value)
    {
        if (value)
        {
            life.boundaryStart = life.boundaryStart+(-1);
            life.boundaryEnd = life.boundaryEnd+(-1);
        }
        else
        {
            life.boundaryStart = life.boundaryStart+(1);
            life.boundaryEnd = life.boundaryEnd+(1);
        }

        return life;
    }
    // 重载加法运算符
    public static Life operator +(Life life1, Life life2)
    {
        long newBoundaryStart = Math.Min(life1.boundaryStart, life2.boundaryStart);
        long newBoundaryEnd = Math.Max(life1.boundaryEnd, life2.boundaryEnd);

        Life newLife = new Life(life1.lifespan, newBoundaryStart, newBoundaryEnd);
        return newLife;
    }
}